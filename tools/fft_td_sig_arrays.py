# -*- coding: utf-8 -*-
##
# \file     fft_td_sig_arrays.py
# \title    Calculation of the FFT for the time domain pressure signals.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 26 Apr.
##
import numpy as np
from scipy import signal
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).
                    split( os.path.sep ) )


def basic_fft(p, t, Ts, time_axis, window):
    """
     Basic FFT that can be used for time-domain pressures.

    :param  p:   acoustic pressure (Pa),  2d or 3d array (dist., angles, time).
    :param  t:   discrete time sequence (s), 1d array.
    :param  Ts:  time step (s), scalar.
    :param  time_axis:  axis of the pressure array that corresponds to
            the time sequence.
    :rtype  time_axis: integer
    :param  window:  name of the window applied to the time signal.
    :rtype  window: string

    :param      Fs  sampling frequency (Hz), scalar.
    :param      n   length of the time sequence, scalar.
    :param      f   frequency sequence (s), 1d array.
    :param      w   hanning window, 1d array.
    :param      P   acoustic pressure as function of frequency, 1d array.

    :return     pressure as a function of the frequency.
    """
    # =========================================================================
    #   Set the parameters for the FFT
    # =========================================================================
    Fs      = 1./Ts
    n       = len(t)
    f       = Fs / 1.0 * np.arange(0, n) / n

    if window == 'rectangular' or window == 'boxcar':
        w = signal.boxcar(p.shape[time_axis])
    elif window == 'hanning':
        w = np.hanning(p.shape[time_axis])
    elif window == 'hamming':
        w = np.hamming(p.shape[time_axis])
    else:
        w = 1.
        print 'No time window applied before FFT.'

    pw = p*w

    # P = np.fft.fft(pw, n=p.shape[time_axis], axis=time_axis)*2. / \
    #     np.sqrt(p.shape[time_axis])
    P = np.fft.fft(pw, n=p.shape[time_axis], axis=time_axis)*1. / \
        p.shape[time_axis]
    # not returning the DC value [1:] to avoid the multiply by zero in miki
    return np.abs(P[1:]), f[1:]


def basic_fft_inner_loading(num_method, case, h_num, time_axis, window):
    """
    FFT for case 1 and the calculation of orders of acc errors_calc3_*
    :param num_method: name of the numerical method of interest.
    :type num_method: string
    :param case: number of the case
    :type case: int
    :param h_num: number of the grid size.
    :type h_num: int
    :param time_axis: axis that corresponds to the time sequence.
    :type time_axis: int
    :param window: name of the window applied to the time signal.
    :type window: string
    :return: the frequency and pressure as a function of frequency
    :rtype: list and arrays of floats
    """
    # =========================================================================
    #   Load the saved pressures and additional variables
    # =========================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, num_method)
    t = np.load(os.path.join(res_path, 't_%i.npy' % h_num))
    Ts = np.load(os.path.join(res_path, 'Ts_%i.npy' % h_num))
    if case == 1:
        p_axial = np.load(os.path.join(res_path, 'p_axial_%i.npy' % h_num))
        p_diag = np.load(os.path.join(res_path, 'p_diag_%i.npy' % h_num))

    # =========================================================================
    #   Set the parameters for the FFT
    # =========================================================================
    Fs = 1./Ts
    n = len(t)
    f = Fs / 1.0 * np.arange(0, n) / n

    if window == 'rectangular' or window == 'boxcar':
        w = signal.boxcar(p.shape[time_axis])
    elif window == 'hanning':
        w = np.hanning(p.shape[time_axis])
    elif window == 'hamming':
        w = np.hamming(p.shape[time_axis])
    else:
        w = 1.
        print 'No time window applied before FFT.'

    if case == 1:
        p_axialw = p_axial * w
        p_diagw = p_diag * w
        P_axialw = np.fft.fft(p_axialw, n=p_axialw.shape[time_axis],
                              axis=time_axis) * 1. / \
                   p_axialw.shape[time_axis]
        P_diagw = np.fft.fft(p_diagw, n=p_diagw.shape[time_axis],
                              axis=time_axis) * 1. / \
                  p_diagw.shape[time_axis]

    # not returning the DC value [1:] to avoid the multiply by zero in miki
    return np.abs(P_axialw[1:]), np.abs(P_diagw[1:]), f[1:]


def fft_conv_arrays(num_method, case, h_num):
    """
     FFTs that return free fiel and total pressure for the propagation above
     a ground.

    :param  num_method   name of the numerical method, string.
    :param  h_num   number of the grid size, scalar.
    :param      t   discrete time sequence (s), 1d array.
    :param      Ts  time step (s), scalar.
    :param      p_t total acoustic pressure (Pa), 1d array.
    :param      p_f free-field acoustic pressure (Pa), 1d array.
    :param      Fs  sampling frequency (Hz), scalar.
    :param      n   length of the time sequence, scalar.
    :param      f   frequency sequence (s), 1d array.
    :param      w   hanning window, 1d array.
    :param      P_T total acoustic pressure as function of frequency, 1d array.
    :param      P_F free-field acoustic pressure as function of freq., 1d array.

    :return     the pressures as a function of the frequency.
    """
    # =========================================================================
    #   Load the saved pressures and additional variables
    # =========================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results','case%i' %case,num_method )
    t       = np.load(os.path.join(res_path,'t_%i.npy' %(h_num)))
    Ts      = np.load(os.path.join(res_path,'Ts_%i.npy' %(h_num)))
    p_t     = np.load(os.path.join(res_path,'p_t_%i.npy' %(h_num)))
    p_f     = np.load(os.path.join(res_path,'p_f_%i.npy' %(h_num)))

    # =========================================================================
    #   Set the parameters for the FFT
    # =========================================================================
    Fs      = 1./Ts
    t       = t[:-1]
    n       = len(t)
    f       = Fs/1.0*np.arange(0, n)/n
    # w       = 1.
    # print 'window: none'
    print np.shape(p_t)
    w       = signal.boxcar(p_t.shape[2])
    # w       = np.hanning(p_t.shape[2])
    ptw     = p_t*w
    pfw     = p_f*w
    # Normalization for one-sided FFT (Parseval th.) and 1/sqrt(N)
    # P_T     = np.fft.fft(ptw, n=p_t.shape[2], axis=2)*2./np.sqrt(p_t.shape[2])
    # P_F     = np.fft.fft(pfw, n=p_f.shape[2], axis=2)*2./np.sqrt(p_f.shape[2])
    # Basic normalization using 1/N
    P_T     = np.fft.fft(ptw, n=p_t.shape[2], axis=2) * 1. / p_t.shape[2]
    P_F     = np.fft.fft(pfw, n=p_f.shape[2], axis=2) * 1. / p_f.shape[2]

    # not returning the DC value [1:] to avoid the multiply by zero in miki
    return np.abs(P_T[:,:,1:]), np.abs(P_F[:,:,1:]), f[1:]


def amplitude_spectrum(x, fs, N=None):
    """
    Basic FFT from the package python-acoustics.
    Amplitude spectrum of instantaneous signal :math:`x(t)`.

    :param x: Instantaneous signal :math:`x(t)`.
    :param fs: Sample frequency :math:`f_s`.
    :param N: Amount of FFT bins.

    The amplitude spectrum gives the amplitudes of the sinusoidal the signal is
    built up from, and the RMS (root-mean-square) amplitudes can easily be found
    by dividing these amplitudes with :math:`\\sqrt{2}`.

    The amplitude spectrum is double-sided.

    """
    N = N if N else x.shape[-1]
    fr = np.fft.fft(x, n=N) / N
    f = np.fft.fftfreq(N, 1.0/fs)
    return np.fft.fftshift(f), np.fft.fftshift(fr, axes=[-1])


def fft_src_calib(num_method, case):
    """
     Basic FFT that can be used for time-domain pressures.

    :param  num_method   name of the numerical method, string.
    :param      t   discrete time sequence (s), 1d array.
    :param      Ts  time step (s), scalar.
    :param      p_t total acoustic pressure (Pa), 1d array.
    :param      p_f free-field acoustic pressure (Pa), 1d array.
    :param      Fs  sampling frequency (Hz), scalar.
    :param      n   length of the time sequence, scalar.
    :param      f   frequency sequence (s), 1d array.
    :param      w   hanning window, 1d array.
    :param      P_T total acoustic pressure as function of frequency, 1d array.
    :param      P_F free-field acoustic pressure as function of frequency, 1d array.

    :return     Return the pressures as a function of the frequency.
    """
    # =========================================================================
    #   Load the saved pressures and additional variables
    # =========================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, num_method)
    t = np.load(os.path.join(res_path, 't.npy'))
    Ts = np.load(os.path.join(res_path, 'Ts.npy'))
    p_f = np.load(os.path.join(res_path, 'p_f.npy'))

    # =========================================================================
    #   Set the parameters for the FFT
    # =========================================================================
    Fs = 1. / Ts
    t = t[:-1]
    n = len(t)
    f = Fs / 1.0 * np.arange(0, n) / n
    w = np.hanning(len(p_f))
    pfw = p_f * w

    P_F = np.fft.fft(pfw, n=len(p_f), axis=0) * 2. / np.sqrt(len(p_f))

    # not returning the DC value [1:] to avoid the multiply by zero in miki
    return P_F[5:], f[5:]