/*
 * ppg_pro.h
 *
 *  Created on: Aug 15, 2025
 *      Author: mathe
 */

#ifndef PPG_PROCESSING_H
#define PPG_PROCESSING_H

#include <stdint.h>

// Filtering functions
void butterworthbpf(int n, double s, double f1, double f2, float frame[], int size);
void butterworthlpf(int n, float s, float f, float waveform[], int size);

// Cycle detection
typedef struct {
    uint32_t last_peak_idx;
    float thresh;
    uint8_t armed;
} CycleDetector;

void CycleDetector_Init(CycleDetector *cd);
int CycleDetector_Process(CycleDetector *cd, float sample, uint32_t idx, float fs, float *ibi_ms);

// Processing pipeline
void PPG_ProcessFrame(float *frame, int size, float fs);

#endif
