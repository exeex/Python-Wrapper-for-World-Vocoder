#include "world/synthesis.h"
#include "synthesis_pulse.h"
#include <math.h>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

static void GetTemporalParametersForTimeBase(const double *f0, int f0_length,
    int fs, int y_length, double frame_period, double lowest_f0,
    double *time_axis, double *coarse_time_axis, double *coarse_f0,
    double *coarse_vuv) {
  for (int i = 0; i < y_length; ++i)
    time_axis[i] = i / static_cast<double>(fs);
  // the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
  for (int i = 0; i < f0_length; ++i) {
    coarse_time_axis[i] = i * frame_period;
    coarse_f0[i] = f0[i] < lowest_f0 ? 0.0 : f0[i];
    coarse_vuv[i] = coarse_f0[i] == 0.0 ? 0.0 : 1.0;
  }
  coarse_time_axis[f0_length] = f0_length * frame_period;
  coarse_f0[f0_length] = coarse_f0[f0_length - 1] * 2 -
    coarse_f0[f0_length - 2];
  coarse_vuv[f0_length] = coarse_vuv[f0_length - 1] * 2 -
    coarse_vuv[f0_length - 2];
}
static int GetPulseLocationsForTimeBase(const double *interpolated_f0,
    const double *time_axis, int y_length, int fs, double *pulse_locations,
    int *pulse_locations_index, double *pulse_locations_time_shift) {
  double *total_phase = new double[y_length];
  double *wrap_phase = new double[y_length];
  double *wrap_phase_abs = new double[y_length - 1];
  total_phase[0] = 2.0 * world::kPi * interpolated_f0[0] / fs;
  wrap_phase[0] = fmod(total_phase[0], 2.0 * world::kPi);
  for (int i = 1; i < y_length; ++i) {
    total_phase[i] = total_phase[i - 1] +
      2.0 * world::kPi * interpolated_f0[i] / fs;
    wrap_phase[i] = fmod(total_phase[i], 2.0 * world::kPi);
    wrap_phase_abs[i - 1] = fabs(wrap_phase[i] - wrap_phase[i - 1]);
  }

  int number_of_pulses = 0;
  for (int i = 0; i < y_length - 1; ++i) {
    if (wrap_phase_abs[i] > world::kPi) {
      pulse_locations[number_of_pulses] = time_axis[i];
      pulse_locations_index[number_of_pulses] = i;

      // calculate the time shift in seconds between exact fractional pulse
      // position and the integer pulse position (sample i)
      // as we don't have access to the exact pulse position, we infer it
      // from the point between sample i and sample i + 1 where the
      // accummulated phase cross a multiple of 2pi
      // this point is found by solving y1 + x * (y2 - y1) = 0 for x, where y1
      // and y2 are the phases corresponding to sample i and i + 1, offset so
      // they cross zero; x >= 0
      double y1 = wrap_phase[i] - 2.0 * world::kPi;
      double y2 = wrap_phase[i + 1];
      double x = -y1 / (y2 - y1);
      pulse_locations_time_shift[number_of_pulses] = x / fs;

      ++number_of_pulses;
    }
  }

  delete[] wrap_phase_abs;
  delete[] wrap_phase;
  delete[] total_phase;

  return number_of_pulses;
}
static int GetTimeBase(const double *f0, int f0_length, int fs,
    double frame_period, int y_length, double lowest_f0,
    double *pulse_locations, int *pulse_locations_index,
    double *pulse_locations_time_shift, double *interpolated_vuv) {
  double *time_axis = new double[y_length];
  double *coarse_time_axis = new double[f0_length + 1];
  double *coarse_f0 = new double[f0_length + 1];
  double *coarse_vuv = new double[f0_length + 1];
  GetTemporalParametersForTimeBase(f0, f0_length, fs, y_length, frame_period,
      lowest_f0, time_axis, coarse_time_axis, coarse_f0, coarse_vuv);
  double *interpolated_f0 = new double[y_length];
  interp1(coarse_time_axis, coarse_f0, f0_length + 1,
      time_axis, y_length, interpolated_f0);
  interp1(coarse_time_axis, coarse_vuv, f0_length + 1,
      time_axis, y_length, interpolated_vuv);

  for (int i = 0; i < y_length; ++i) {
    interpolated_vuv[i] = interpolated_vuv[i] > 0.5 ? 1.0 : 0.0;
    interpolated_f0[i] =
      interpolated_vuv[i] == 0.0 ? world::kDefaultF0 : interpolated_f0[i];
  }

  int number_of_pulses = GetPulseLocationsForTimeBase(interpolated_f0,
      time_axis, y_length, fs, pulse_locations, pulse_locations_index,
      pulse_locations_time_shift);

  delete[] coarse_vuv;
  delete[] coarse_f0;
  delete[] coarse_time_axis;
  delete[] time_axis;
  delete[] interpolated_f0;

  return number_of_pulses;
}


void Synthesis_pulse(const double *f0, int f0_length,
//    const double * const *spectrogram, const double * const *aperiodicity,
    int fft_size, double frame_period, int fs, int y_length, double *y) {


  for (int i = 0; i < y_length; ++i) y[i] = 0.0;

  double *pulse_locations = new double[y_length];
  int *pulse_locations_index = new int[y_length];
  double *pulse_locations_time_shift = new double[y_length];
  double *interpolated_vuv = new double[y_length];
  int number_of_pulses = GetTimeBase(f0, f0_length, fs, frame_period / 1000.0,
      y_length, fs / fft_size + 1.0, pulse_locations, pulse_locations_index,
      pulse_locations_time_shift, interpolated_vuv);


  frame_period /= 1000.0;
  int noise_size;
  int index, offset, lower_limit, upper_limit;
  for (int i = 0; i < number_of_pulses; ++i) {

    offset = pulse_locations_index[i] - fft_size / 2 + 1;

    if(0<=offset && offset < y_length)    y[offset] = 1.0;
  }

  delete[] pulse_locations;
  delete[] pulse_locations_index;
  delete[] pulse_locations_time_shift;
  delete[] interpolated_vuv;
}

void Synthesis_pulse_new(const double *f0, int f0_length, int f_n,
    int fft_size, double frame_period, int fs, int y_length, int *y) {

  // multiple f0 by f_n

  double * fz = new double[f0_length];

  for (int i = 0; i < f0_length; ++i) fz[i] = f0[i]*f_n;

  // init y
  for (int i = 0; i < y_length; ++i) y[i] = 0;

  // make fn one-hot
  if(f_n < 0 || f_n > 32) return; // check 0 < f_n < 32
  int f_n_one_hot[32];
  for(int i=0; i<32; i++) f_n_one_hot[i] = 1 << i;


  double *pulse_locations = new double[y_length];
  int *pulse_locations_index = new int[y_length];
  double *pulse_locations_time_shift = new double[y_length];
  double *interpolated_vuv = new double[y_length];
  int number_of_pulses = GetTimeBase(fz, f0_length, fs, frame_period / 1000.0,
      y_length, fs / fft_size + 1.0, pulse_locations, pulse_locations_index,
      pulse_locations_time_shift, interpolated_vuv);


  frame_period /= 1000.0;
  int noise_size;
  int index, offset, lower_limit, upper_limit;
  for (int i = 0; i < number_of_pulses; ++i) {

    offset = pulse_locations_index[i] - fft_size / 2 + 1;

    if(0<=offset && offset < y_length)    y[offset] = f_n_one_hot[i % f_n];
  }

  delete[] pulse_locations;
  delete[] pulse_locations_index;
  delete[] pulse_locations_time_shift;
  delete[] interpolated_vuv;
}