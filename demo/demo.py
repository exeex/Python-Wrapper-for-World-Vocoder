from __future__ import division, print_function

import os
from shutil import rmtree
import argparse

import numpy as np

# import matplotlib      # Remove this line if you don't need them
# matplotlib.use('Agg')  # Remove this line if you don't need them
import matplotlib.pyplot as plt

import soundfile as sf
# import librosa
import pyworld_mod.pyworld as pw

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--frame_period", type=float, default=5.0)
parser.add_argument("-s", "--speed", type=int, default=1)

EPSILON = 1e-8


def main():
    x, fs = sf.read('utterance/vaiueo2d.wav')

    # 2-3 Harvest with F0 refinement (using Stonemask)

    f0, t = pw.harvest(x, fs)
    f0 = pw.stonemask(x, f0, t, fs)
    sp = pw.cheaptrick(x, f0, t, fs)
    ap = pw.d4c(x, f0, t, fs)
    # y = pw.synthesize(f0, sp, ap, fs, pw.default_frame_period)
    # y_pulse = pw.synthesize_pulse(f0, sp, ap, fs, pw.default_frame_period)

    # y_pulse = pw.synthesize_pulse(f0, fs, frame_period=pw.default_frame_period)
    y_pulse = pw.synthesize_pulse(pw.fill_0hz(f0), fs, fft_size=1024, frame_period=pw.default_frame_period)
    print(sp.shape)
    print('Please check "test" directory for output files')
    return y_pulse, f0

# def fill_0hz(f0):
#     zeros = np.where(f0 == 0.0)[0]
#     if zeros[1] == 0:
#         zeros = zeros[1:]
#     for zero_idx in zeros:
#         f0[zero_idx] = f0[zero_idx - 1]
#     return f0

if __name__ == '__main__':
    args = parser.parse_args()
    x, fs = sf.read('utterance/vaiueo2d.wav')
    y_pulse, f0 = main()

    # zeros = np.where(f0 == 0.0)[0]
    # sample_per_frame = fs / (1000 / pw.default_frame_period)
    #
    # for frame_idx in zeros:
    #     sample_start_idx = int(frame_idx * sample_per_frame)
    #     sample_end_idx = int((frame_idx + 1) * sample_per_frame)
    #     y_pulse[sample_start_idx:sample_end_idx] = 0




    plt.plot(f0)
    plt.show()
    plt.plot(y_pulse)
    plt.show()


