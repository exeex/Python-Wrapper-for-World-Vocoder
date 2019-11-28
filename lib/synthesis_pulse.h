#ifndef WORLD_SYNTHESIS_PULSE_H_
#define WORLD_SYNTHESIS_PULSE_H_

#include "world/macrodefinitions.h"


void Synthesis_pulse(const double *f0, int f0_length,
//    const double * const *spectrogram, const double * const *aperiodicity,
    int fft_size, double frame_period, int fs, int y_length, double *y);


#endif