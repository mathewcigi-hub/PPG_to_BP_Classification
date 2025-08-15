/*
 * PPG_features.h
 *
 *  Created on: Aug 29, 2024
 *      Author: mathe
 */

#ifndef PPG_FEATURES_H
#define PPG_FEATURES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// Public constants
#define PPGF_OK                 0
#define PPGF_ERR_BOUNDS       (-1)
#define PPGF_ERR_SHORT        (-2)

// Feature vector (36 total)
// - First 31 are “values” (slopes, areas, intensities, times, ratios)
// - Last 5 are “locations” (sample indices) useful for QC / debugging
typedef struct {
    // 1..31: value features
    float Sys_time;     // Peak_Systolic_L - ED1_L
    float Sys_I_diff;   // Pulse[PS] - Pulse[ED1]
    float D_AS;         // (Pulse[PS]-Pulse[ED1])/(PS-ED1)
    float Di_time;      // ED2 - PS
    float Di_I_diff;    // Pulse[ED2] - Pulse[PS]
    float D_DS;         // Di_I_diff / Di_time
    float FD_Rise_T;    // FDM - FDV1
    float FD_Rise_I;    // FD[FDM] - FD[FDV1]
    float FD_AS;        // FD_Rise_I / FD_Rise_T
    float FD_Dec_T;     // FDV2 - FDM
    float FD_Dec_I;     // FD[FDV2] - FD[FDM]
    float FD_DS;        // FD_Dec_I / FD_Dec_T
    float SD_Rise_T;    // SDM - SDV1
    float SD_Rise_I;    // SD[SDM] - SD[SDV1]
    float SD_AS;        // SD_Rise_I / SD_Rise_T
    float SD_Dec_T;     // SDV2 - SDM
    float SD_Dec_I;     // SD[SDV2] - SD[SDM]
    float SD_DS;        // SD_Dec_I / SD_Dec_T
    float FD_AA;        // sum(FD - FD[FDV1], FDV1..FDM)
    float FD_DA;        // sum(FD - FD[FDV2], FDM..FDV2)
    float SD_AA;        // sum(SD - SD[SDV1], SDV1..SDM)
    float SD_DA;        // sum(SD - SD[SDV2], SDM..SDV2)
    float SO_Pk_S;      // (FoD[SO] - Pulse[PS])/(SO - PS)
    float DN_Pk_S;      // (Pulse[DN] - Pulse[PS])/(DN - PS)
    float Di_D;         // Pulse[ED1]
    float Max_dD;       // (Baseline-corrected delta) at PS
    float dD_T;         // Max_dD / Di_D
    float SO_Pk_T;      // SO - PS
    float DN_Mi_T;      // ED2 - DN
    float T;            // ED2 - ED1
    float Fmin;         // ED_Amp (= Pulse[ED1] in this impl)

    // 32..36: locations (indices) — integer-like, but stored as float
    float L_ED1;
    float L_ED2;
    float L_PS;     // Peak_Systolic_L
    float L_DN;     // diastolic notch
    float L_SO;     // zero-cross based point on 4th deriv
} PPGFeatures;

// Main API
// Inputs:
//   pulse[] : single PPG cycle (float), length N (>= 20 recommended)
//   N       : number of samples in the cycle
//   fs      : sampling rate (Hz), e.g., 125
// Outputs:
//   out     : filled feature vector (36 fields)
// Returns PPGF_OK on success, negative on error.
int PPG_ExtractFeatures(const float *pulse, int N, float fs, PPGFeatures *out);

#ifdef __cplusplus
}
#endif

#endif // PPG_FEATURES_H
