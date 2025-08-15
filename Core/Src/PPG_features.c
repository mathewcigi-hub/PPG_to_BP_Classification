#include "PPG_features.h"
#include <string.h>   // memset
#include <math.h>     // fabsf, fmaxf, fminf

// ---------- helpers ----------
static inline int clampi(int x, int lo, int hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}
static inline float fsafe_div(float num, float den) {
    return (fabsf(den) > 1e-6f) ? (num / den) : 0.0f;
}

// simple first/second/fourth discrete derivatives (no noise shaping)
static void diff1(const float *x, int N, float *dx) {
    if (N < 2) return;
    dx[0] = 0.0f;
    for (int i = 1; i < N; ++i) dx[i] = x[i] - x[i-1];
}
static void diff2(const float *x, int N, float *d2x) {
    if (N < 3) { for (int i=0;i<N;i++) d2x[i]=0; return; }
    d2x[0] = 0.0f;
    for (int i = 1; i < N-1; ++i) d2x[i] = x[i+1] - 2.0f*x[i] + x[i-1];
    d2x[N-1] = 0.0f;
}
static void diff4(const float *x, int N, float *d4x) {
    // apply 2nd derivative twice
    static float tmp[4096];
    if (N > (int)(sizeof(tmp)/sizeof(tmp[0]))) {
        // fallback: zero if too long (you can increase size as needed)
        for (int i=0;i<N;i++) d4x[i]=0.0f;
        return;
    }
    diff2(x, N, tmp);
    diff2(tmp, N, d4x);
}

// local extrema: argrelextrema like (strict)
static int find_local_minima(const float *x, int N, int *idx, int max_idx) {
    int k = 0;
    for (int i = 1; i < N-1 && k < max_idx; ++i) {
        if (x[i] < x[i-1] && x[i] < x[i+1]) idx[k++] = i;
    }
    return k;
}
static int find_local_maxima(const float *x, int N, int *idx, int max_idx) {
    int k = 0;
    for (int i = 1; i < N-1 && k < max_idx; ++i) {
        if (x[i] > x[i-1] && x[i] > x[i+1]) idx[k++] = i;
    }
    return k;
}
static int argmax(const float *x, int N) {
    if (N <= 0) return 0;
    int m = 0; float mv = x[0];
    for (int i=1;i<N;i++){ if (x[i] > mv) { mv = x[i]; m = i; } }
    return m;
}
static int argmin(const float *x, int N) {
    if (N <= 0) return 0;
    int m = 0; float mv = x[0];
    for (int i=1;i<N;i++){ if (x[i] < mv) { mv = x[i]; m = i; } }
    return m;
}

// area helper: sum_{i=a..b} (x[i] - base)
static float area_minus_base(const float *x, int a, int b, float base) {
    if (b < a) return 0.0f;
    float s = 0.0f;
    for (int i=a; i<=b; ++i) s += (x[i] - base);
    return s;
}

// baseline correction (simple ramp index from start to end minima, per your sketch)
static void baseline_correct(const float *pulse, int N, float minima_offset,
                             float *dD_out, float *D_out) {
    int half = N/2;
    int end_i   = argmin(pulse, half > 0 ? half : N);
    int start_i = half < N ? argmin(pulse + half, N - half) + half : 0;

    // linear ramp from start->end across N samples (index-space)
    float rampN = (float)(N - 1);
    float slope = (end_i - start_i) / (rampN > 0 ? rampN : 1.0f);

    for (int i = 0; i < N; ++i) {
        float ramp = start_i + slope * i; // "index ramp"
        dD_out[i] = pulse[i] - ramp;
        D_out[i]  = dD_out[i] + minima_offset;
    }
}

// DN & SDMax from SD up to ED2, enforce min separation (fs/10)
static int find_DN_SDMax(const float *SD, int N, int ED2, float fs, int *DN_L, int *SDMax_L) {
    int mins[256]; int nmins = find_local_minima(SD, (ED2>0 && ED2<N) ? ED2 : N, mins, 256);
    int maxs[256]; int nmaxs = find_local_maxima(SD, (ED2>0 && ED2<N) ? ED2 : N, maxs, 256);

    // keep only positive maxima (as in python: array[peaks]>0)
    int pos_maxs[256]; int npos = 0;
    for (int i=0;i<nmaxs && npos<256;i++){
        if (SD[maxs[i]] > 0.0f) pos_maxs[npos++] = maxs[i];
    }
    if (npos < 2) { *DN_L = 0; *SDMax_L = (npos>0?pos_maxs[0]:0); return 0; }

    int min_sep = (int)(fs/10.0f); if (min_sep < 1) min_sep = 1;

    for (int i=0;i<npos-1;i++){
        if (pos_maxs[i+1] - pos_maxs[i] > min_sep) {
            *DN_L = pos_maxs[i+1];
            *SDMax_L = pos_maxs[i];
            return 1;
        }
    }
    // fallback
    *SDMax_L = pos_maxs[0];
    *DN_L    = pos_maxs[npos-1];
    return 0;
}

// zero crossing after SDMax on -SD (first > 0)
static int second_der_zero_cross(const float *SD, int N, int SDMax_L) {
    if (SDMax_L < 0 || SDMax_L >= N) return 0;
    for (int i = SDMax_L; i < N; ++i) {
        if (-SD[i] > 0.0f) return i;
    }
    return SDMax_L;
}

// 4th-deriv “SO” (as in your sketch):
// 1) find first 4th-deriv local max after SDMax -> FoD_peak
// 2) from FoD_peak forward, find first i where -FoD[i] > 0 -> SO
static int fourth_der_so(const float *FoD, int N, int SDMax_L) {
    int peaks[256]; int npeaks = find_local_maxima(FoD, N, peaks, 256);
    int FoD_peak = SDMax_L;
    for (int i=0;i<npeaks;i++){
        if (peaks[i] > SDMax_L) { FoD_peak = peaks[i]; break; }
    }
    for (int i = FoD_peak; i < N; ++i){
        if (-FoD[i] > 0.0f) return i;
    }
    return FoD_peak;
}

// first/second derivative valley logic restricted to first half
static void fd_sd_valleys(const float *D, int N, int *M_L, int *V1_L, int *V2_L) {
    int half = (N>0)? (N/2) : 0;
    if (half < 3) { *M_L = *V1_L = *V2_L = 0; return; }

    int M = argmax(D, half);
    int mins[256]; int nmins = find_local_minima(D, half, mins, 256);

    int v1 = 0, v2 = 0; int found_v1 = 0, found_v2 = 0;
    for (int i=0;i<nmins;i++){
        if (mins[i] > M) { v1 = mins[i]; found_v1 = 1; break; }
    }
    for (int i=0;i<nmins;i++){
        if (mins[i] < M) { v2 = mins[i]; found_v2 = 1; } // last one before M
    }
    *M_L  = M;
    *V1_L = found_v1 ? v1 : 0;
    *V2_L = found_v2 ? v2 : 0;
}

// intersecting tangent (your sketch): IT_point = (Pulse_min1 - c)/MaxFD,
// where c = Pulse_at(MaxFD_idx) - MaxFD_idx*MaxFD_val, using first third.
static float intersect_tangent(const float *FD, const float *Pulse, int N) {
    int third = (N>0)? (N/3) : 0;
    if (third < 3) return 0.0f;
    int m_idx = argmax(FD, third);
    float m_val = FD[m_idx];
    float c = Pulse[m_idx] - (m_idx * m_val);

    int min_idx = argmax(Pulse, third); // NOTE: original code likely intended argmin; kept as given
    float min_val = Pulse[min_idx];
    float S1 = min_val - c;

    return fsafe_div(S1, m_val);
}

// ---------- main API ----------
int PPG_ExtractFeatures(const float *pulse, int N, float fs, PPGFeatures *out) {
    if (!pulse || !out) return PPGF_ERR_BOUNDS;
    if (N < 20) { memset(out, 0, sizeof(*out)); return PPGF_ERR_SHORT; }

    // Work buffers (no malloc) — adjust max length if needed
    if (N > 4096) return PPGF_ERR_BOUNDS;
    static float FD[4096], SD[4096], FoD[4096];
    static float dD[4096], D[4096];

    // Build derivatives
    diff1(pulse, N, FD);
    diff2(pulse, N, SD);
    diff4(pulse, N, FoD);

    // Key points from pulse
    int PS = argmax(pulse, N);

    // valleys across full cycle; choose ED2 = last valley; ED1 = last valley before PS
    int v[512]; int nv = find_local_minima(pulse, N, v, 512);
    int ED1 = 0, ED2 = (nv>0? v[nv-1] : (N-1));
    for (int i=0;i<nv;i++){
        if (v[i] < PS) ED1 = v[i];
        else break;
    }
    if (ED1 >= PS) ED1 = (PS>0? PS-1 : 0); // sanity

    // delD & ED_Amp per your notes (ED_Amp = pulse[ED1])
    int half = N/2;
    int fh_N = (half>0? half : N);
    int fh_max = argmax(pulse, fh_N);
    int fh_min = argmin(pulse, fh_N);
    float delD = pulse[fh_max] - pulse[fh_min];
    float ED_Amp = pulse[ED1];

    // First/Second derivative landmarks in first half
    int FDM=0, FDV1=0, FDV2=0; fd_sd_valleys(FD, N, &FDM, &FDV1, &FDV2);
    int SDM=0, SDV1=0, SDV2=0; fd_sd_valleys(SD, N, &SDM, &SDV1, &SDV2);

    // SD: DN and SDMax up to ED2
    int DN=0, SDMax=0; (void)find_DN_SDMax(SD, N, ED2, fs, &DN, &SDMax);

    // 2nd-der zero cross (flow_peak) and 4th-der “SO”
    int Flow_peak = second_der_zero_cross(SD, N, SDMax);
    int SO = fourth_der_so(FoD, N, SDMax);

    // Intersecting tangent (kept as in your sketch)
    float IT_point = intersect_tangent(FD, pulse, N);
    (void)IT_point; // not part of the 36 exported, but here if you want to add

    // Baseline correction
    baseline_correct(pulse, N, 0.0f /*Minima*/, dD, D);

    // ===== Features =====
    PPGFeatures f; memset(&f, 0, sizeof(f));

    f.Sys_time   = (float)(PS - ED1);
    f.Sys_I_diff = pulse[PS] - pulse[ED1];
    f.D_AS       = fsafe_div(f.Sys_I_diff, f.Sys_time);

    f.Di_time    = (float)(ED2 - PS);
    f.Di_I_diff  = pulse[ED2] - pulse[PS];
    f.D_DS       = fsafe_div(f.Di_I_diff, f.Di_time);

    f.FD_Rise_T  = (float)(FDM - FDV1);
    f.FD_Rise_I  = FD[FDM] - FD[FDV1];
    f.FD_AS      = fsafe_div(f.FD_Rise_I, f.FD_Rise_T);

    f.FD_Dec_T   = (float)(FDV2 - FDM);
    f.FD_Dec_I   = FD[FDV2] - FD[FDM];
    f.FD_DS      = fsafe_div(f.FD_Dec_I, f.FD_Dec_T);

    f.SD_Rise_T  = (float)(SDM - SDV1);
    f.SD_Rise_I  = SD[SDM] - SD[SDV1];
    f.SD_AS      = fsafe_div(f.SD_Rise_I, f.SD_Rise_T);

    f.SD_Dec_T   = (float)(SDV2 - SDM);
    f.SD_Dec_I   = SD[SDV2] - SD[SDM];
    f.SD_DS      = fsafe_div(f.SD_Dec_I, f.SD_Dec_T);

    // Areas
    f.FD_AA      = area_minus_base(FD, clampi(FDV1,0,N-1), clampi(FDM,0,N-1), FD[FDV1]);
    f.FD_DA      = area_minus_base(FD, clampi(FDM,0,N-1), clampi(FDV2,0,N-1), FD[FDV2]);
    f.SD_AA      = area_minus_base(SD, clampi(SDV1,0,N-1), clampi(SDM,0,N-1), SD[SDV1]);
    f.SD_DA      = area_minus_base(SD, clampi(SDM,0,N-1), clampi(SDV2,0,N-1), SD[SDV2]);

    // Slopes vs PS
    f.SO_Pk_S    = fsafe_div(FoD[SO] - pulse[PS], (float)(SO - PS));
    f.DN_Pk_S    = fsafe_div(pulse[DN] - pulse[PS], (float)(DN - PS));

    // Baseline-derived features
    f.Di_D       = pulse[ED1];
    f.Max_dD     = dD[PS];
    f.dD_T       = fsafe_div(f.Max_dD, f.Di_D);

    // Times
    f.SO_Pk_T    = (float)(SO - PS);
    f.DN_Mi_T    = (float)(ED2 - DN);
    f.T          = (float)(ED2 - ED1);

    // “Fmin” = ED_Amp
    f.Fmin       = ED_Amp;

    // Locations
    f.L_ED1 = (float)ED1;
    f.L_ED2 = (float)ED2;
    f.L_PS  = (float)PS;
    f.L_DN  = (float)DN;
    f.L_SO  = (float)SO;

    // Basic sanity: zero any NAN/INF
    const float *pf = (const float*)&f;
    for (int i=0;i<(int)(sizeof(PPGFeatures)/sizeof(float));++i){
        float v = pf[i];
        if (!(v == v) || v > 1e20f || v < -1e20f) {
            ((float*)&f)[i] = 0.0f;
        }
    }

    *out = f;
    (void)delD; // available if you want to add as a feature
    return PPGF_OK;
}
