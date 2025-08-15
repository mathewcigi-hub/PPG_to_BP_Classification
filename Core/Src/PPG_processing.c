/*
 * PPG_processing.c
 *
 *  Created on: Aug 11, 2024
 *      Author: mathe
 */


#include "ppg_processing.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* =====================================================
   Filtering
   ===================================================== */
void butterworthbpf(int n, double s, double f1, double f2, float frame[], int size) {
    double a = cos(M_PI * (f1 + f2) / s) / cos(M_PI * (f1 - f2) / s);
    double a2 = a * a;
    double b = tan(M_PI * (f1 - f2) / s);
    double b2 = b * b;
    double r;
    int i, j;

    n = n / 4;
    double *A  = malloc(n * sizeof(double));
    double *d1 = malloc(n * sizeof(double));
    double *d2 = malloc(n * sizeof(double));
    double *d3 = malloc(n * sizeof(double));
    double *d4 = malloc(n * sizeof(double));
    double *w0 = calloc(n, sizeof(double));
    double *w1 = calloc(n, sizeof(double));
    double *w2 = calloc(n, sizeof(double));
    double *w3 = calloc(n, sizeof(double));
    double *w4 = calloc(n, sizeof(double));

    for (i = 0; i < n; ++i) {
        r = sin(M_PI * (2.0 * i + 1.0) / (4.0 * n));
        s = b2 + 2.0 * b * r + 1.0;
        A[i]  = b2 / s;
        d1[i] = 4.0 * a * (1.0 + b * r) / s;
        d2[i] = 2.0 * (b2 - 2.0 * a2 - 1.0) / s;
        d3[i] = 4.0 * a * (1.0 - b * r) / s;
        d4[i] = -(b2 - 2.0 * b * r + 1.0) / s;
    }

    // forward filtering
    for (j = 0; j < size; j++) {
        for (i = 0; i < n; ++i) {
            w0[i] = d1[i] * w1[i] + d2[i] * w2[i] + d3[i] * w3[i] + d4[i] * w4[i] + frame[j];
            frame[j] = A[i] * (w0[i] - 2.0 * w2[i] + w4[i]);
            w4[i] = w3[i]; w3[i] = w2[i]; w2[i] = w1[i]; w1[i] = w0[i];
        }
    }

    // reset
    memset(w0, 0, n * sizeof(double));
    memset(w1, 0, n * sizeof(double));
    memset(w2, 0, n * sizeof(double));
    memset(w3, 0, n * sizeof(double));
    memset(w4, 0, n * sizeof(double));

    // reverse filtering
    for (j = size - 1; j >= 0; j--) {
        for (i = 0; i < n; ++i) {
            w0[i] = d1[i] * w1[i] + d2[i] * w2[i] + d3[i] * w3[i] + d4[i] * w4[i] + frame[j];
            frame[j] = A[i] * (w0[i] - 2.0 * w2[i] + w4[i]);
            w4[i] = w3[i]; w3[i] = w2[i]; w2[i] = w1[i]; w1[i] = w0[i];
        }
    }

    free(A); free(d1); free(d2); free(d3); free(d4);
    free(w0); free(w1); free(w2); free(w3); free(w4);
}

void butterworthlpf(int n, float s, float f, float waveform[], int size) {
    int i, j;
    n = n / 2;
    float a = tan(M_PI * f / s);
    float a2 = a * a;
    float r;
    float *A  = malloc(n * sizeof(float));
    float *d1 = malloc(n * sizeof(float));
    float *d2 = malloc(n * sizeof(float));
    float *w0 = calloc(n, sizeof(float));
    float *w1 = calloc(n, sizeof(float));
    float *w2 = calloc(n, sizeof(float));

    for (i = 0; i < n; ++i) {
        r = sin(M_PI * (2.0 * i + 1.0) / (4.0 * n));
        s = a2 + 2.0 * a * r + 1.0;
        A[i]  = a2 / s;
        d1[i] = 2.0 * (1 - a2) / s;
        d2[i] = -(a2 - 2.0 * a * r + 1.0) / s;
    }

    for (j = 0; j < size; j++) {
        for (i = 0; i < n; ++i) {
            w0[i] = d1[i] * w1[i] + d2[i] * w2[i] + waveform[j];
            waveform[j] = A[i] * (w0[i] + 2.0 * w1[i] + w2[i]);
            w2[i] = w1[i]; w1[i] = w0[i];
        }
    }
    for (j = size - 1; j >= 0; j--) {
        for (i = 0; i < n; ++i) {
            w0[i] = d1[i] * w1[i] + d2[i] * w2[i] + waveform[j];
            waveform[j] = A[i] * (w0[i] + 2.0 * w1[i] + w2[i]);
            w2[i] = w1[i]; w1[i] = w0[i];
        }
    }

    free(A); free(d1); free(d2); free(w0); free(w1); free(w2);
}

/* =====================================================
   Cycle Detection
   ===================================================== */
void CycleDetector_Init(CycleDetector *cd) {
    cd->last_peak_idx = 0;
    cd->thresh = 100.0f;
    cd->armed = 0;
}

int CycleDetector_Process(CycleDetector *cd, float sample, uint32_t idx, float fs, float *ibi_ms) {
    static float prev = 0;
    static float prev_d = 0;
    float d = sample - prev;

    uint32_t min_dist = (uint32_t)(0.3f * fs); // 300 ms refractory
    int found = 0;

    if (!cd->armed && (sample > cd->thresh) && d > 0 && (idx - cd->last_peak_idx) > min_dist) {
        cd->armed = 1;
    }

    if (cd->armed && prev_d > 0 && d <= 0) {
        found = 1;
        if (cd->last_peak_idx > 0) {
            uint32_t ibi_samples = idx - cd->last_peak_idx;
            *ibi_ms = (ibi_samples * 1000.0f) / fs;
        }
        cd->last_peak_idx = idx;
        cd->thresh = 0.8f * cd->thresh + 0.2f * fabsf(sample);
        cd->armed = 0;
    }

    prev_d = d;
    prev = sample;
    return found;
}

/* =====================================================
   Processing Pipeline
   ===================================================== */
void PPG_ProcessFrame(float *frame, int size, float fs) {
    // Example: Band-pass 0.5 - 8 Hz
    butterworthbpf(4, fs, 0.5, 8.0, frame, size);
}
