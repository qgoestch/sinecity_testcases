# -*- coding: utf-8 -*-
##
# \file     level_calc_paumond_facade.py
# \title    Calculation of the levels
# \author   Pierre Chobeau
# \version  0.2
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2018, 10 Jan.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from fft_td_sig_arrays import fft_conv_arrays, basic_fft
from acoustics.signal import third_octaves
import matplotlib.pyplot as plt


def levels(d_sr, d_sr_ref, case, method):
    """
    Calculation of the levels for the receivers defined on the segment d_sr.

    :param d_sr: horizontal distances between the facade and the receivers (m).
    :type d_sr: list of floats
    :param d_sr_ref: horizontal distances between the facade and the pedestrian (m).
    :type d_sr_ref: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param method: name of the numerical method used for the pressure calculation.
    :type method: string
    :return: the levels in dB SPL
    :rtype: list of floats (maybe specrums ?)
    """
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, method)
    t  = np.load(os.path.join(res_path, 't.npy'))
    Ts = np.load(os.path.join(res_path, 'Ts.npy'))
    p0 = np.load(os.path.join(res_path, 'p0.npy'))
    p1 = np.load(os.path.join(res_path, 'p1.npy'))
    p2 = np.load(os.path.join(res_path, 'p2.npy'))
    p3 = np.load(os.path.join(res_path, 'p3.npy'))

    P0, f = zip(*(basic_fft(p0[i, :], t, Ts, 0, 'hanning')
                  for i in range(len(d_sr_ref))))
    P0 = np.asarray(P0)
    l_db_spl1 = [20. * np.log10(np.sum(np.abs(P0[i, :])) / (2. * 10**-5))
                for i in range(len(d_sr_ref))]

    P1, f = zip(*(basic_fft(p1[i, :], t, Ts, 0, 'hanning')
                  for i in range(len(d_sr))))
    P1 = np.asarray(P1)
    l_db_spl1 = [20. * np.log10(np.sum(np.abs(P1[i, :])) / (2. * 10**-5))
                for i in range(len(d_sr))]

    P2, f = zip(*(basic_fft(p2[i, :], t, Ts, 0, 'hanning')
                  for i in range(len(d_sr))))
    P2 = np.asarray(P2)
    l_db_spl2 = [20. * np.log10(np.sum(np.abs(P2[i, :])) / (2. * 10**-5))
                for i in range(len(d_sr))]

    P3, f = zip(*(basic_fft(p3[i, :], t, Ts, 0, 'hanning')
                  for i in range(len(d_sr))))
    P3 = np.asarray(P3)
    l_db_spl3 = [20. * np.log10(np.sum(np.abs(P3[i, :])) / (2. * 10**-5))
                for i in range(len(d_sr))]


    l_db_spl_f_0 = np.zeros((len(d_sr), P0.shape[1]), dtype=np.float64)
    l_db_spl_f_1 = np.zeros((len(d_sr), P0.shape[1]), dtype=np.float64)
    l_db_spl_f_2 = np.zeros((len(d_sr), P0.shape[1]), dtype=np.float64)
    l_db_spl_f_3 = np.zeros((len(d_sr), P0.shape[1]), dtype=np.float64)
    for i in range(len(d_sr)):
        l_db_spl_f_0[i, :] = 20. * np.log10(np.abs(P0[i, :]) / (2. * 10 ** -5))
        l_db_spl_f_1[i, :] = 20. * np.log10(np.abs(P1[i, :]) / (2. * 10 ** -5))
        l_db_spl_f_2[i, :] = 20. * np.log10(np.abs(P2[i, :]) / (2. * 10 ** -5))
        l_db_spl_f_3[i, :] = 20. * np.log10(np.abs(P3[i, :]) / (2. * 10 ** -5))

    l_db_spl_subtract1 = np.zeros((len(d_sr), len(d_sr), P0.shape[1]),
                                  dtype=np.float64)
    l_db_spl_subtract2 = np.zeros((len(d_sr), len(d_sr), P0.shape[1]),
                                  dtype=np.float64)
    l_db_spl_subtract3 = np.zeros((len(d_sr), len(d_sr), P0.shape[1]),
                                  dtype=np.float64)
    for i in range(len(d_sr)):
        for j in range(len(d_sr)):
            l_db_spl_subtract1[i, j, :] = l_db_spl_f_1[j, :] - \
                                          l_db_spl_f_0[i, :]
            l_db_spl_subtract2[i, j, :] = l_db_spl_f_2[j, :] - \
                                          l_db_spl_f_0[i, :]
            l_db_spl_subtract3[i, j, :] = l_db_spl_f_3[j, :] - \
                                          l_db_spl_f_0[i, :]

    print np.shape(l_db_spl_f_1), np.shape(l_db_spl_f_2), np.shape(f)
    f = np.asarray(f)
    # print f

    frequencies = [63., 80., 100., 125., 160., 200.,
                   250., 315., 400., 500., 630., 800., 1000., 1250., 1600,2000., 2500]
    level0_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level1_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level2_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level3_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level1_0_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level2_0_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level3_0_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    p_0_1 = np.zeros((len(d_sr), l_db_spl_subtract1.shape[2]), dtype=np.float64)
    for i in range(len(d_sr)):
        f_third, level0_third[i, :] = third_octaves(p0[i, :], 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        f_third, level1_third[i, :] = third_octaves(p1[i, :], 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        f_third, level2_third[i, :] = third_octaves(p2[i, :], 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        f_third, level3_third[i, :] = third_octaves(p3[i, :], 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        # print np.shape(10 ** (l_db_spl_subtract1[i, 0, :] / 20.)), np.shape(p_0_1[i, :])
        # p_0_1[i, :] = 10 ** (l_db_spl_subtract1[i, 0, :] / 20.)

        f_third, level1_0_third[i, :] = third_octaves(10 ** (l_db_spl_subtract2[i, 0, :]/20.), 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        f_third, level2_0_third[i, :] = third_octaves(10 ** (l_db_spl_subtract2[i, 0, :]/20.), 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        f_third, level3_0_third[i, :] = third_octaves(10 ** (l_db_spl_subtract3[i, 0, :]/20.), 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
    frequencies = np.asarray(frequencies)
    fig = plt.figure('1')
    ax = fig.add_subplot(311)
    plt.semilogx(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies,
           np.median(level2_third[0, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies + 1.5,
           np.median(level3_third[0, :] - level0_third[:, :], axis=0))
    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    plt.legend(('h_r = 4 m', 'h_r = 7 m', 'h_r = 12 m'))
    # plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.grid()

    ax = fig.add_subplot(312)
    plt.semilogx(frequencies-1.5, np.median(level1_third[1, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies,
           np.median(level2_third[1, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies + 1.5,
           np.median(level3_third[1, :] - level0_third[:, :], axis=0))
    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    # plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.grid()

    ax = fig.add_subplot(313)
    plt.semilogx(frequencies-1.5, np.median(level1_third[2, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies,
           np.median(level2_third[2, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies + 1.5,
           np.median(level3_third[2, :] - level0_third[:, :], axis=0))
    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.grid()

    plt.figure('Mic pieton mobile')
    # avg_sub = np.median(np.median(l_db_spl_subtract1, axis=0), axis=0)
    for i in range(len(d_sr)):
        plt.plot(f[0, :], l_db_spl_subtract2[i, 0, :])
    plt.xlabel('f (Hz)', fontsize=16)
    plt.ylabel('Pressure level dB SPL', fontsize=14)
    plt.xlim(0, 1100)
    plt.ylim(-40, 40)
    plt.grid()

    fig = plt.figure('Mic. facade mobile')
    avg_sub = np.median(np.median(l_db_spl_subtract1, axis=0), axis=0)
    ax = fig.add_subplot(311)
    for j in range(len(d_sr)):
        plt.plot(f[0, :], l_db_spl_subtract2[2, j, :])
    # plt.plot(f[0, :], avg_sub[:], 'ro')
    # plt.xlabel('f (Hz)', fontsize=16)
    plt.ylabel('Pressure level dB SPL', fontsize=14)
    plt.title('1er etage : h_r = 2.50m')
    plt.xlim(0, 1200)
    plt.ylim(-40, 40)
    plt.grid()

    ax = fig.add_subplot(312)
    for j in range(len(d_sr)):
        plt.plot(f[0, :], l_db_spl_subtract2[2, j, :])
    # plt.plot(f[0, :], avg_sub[:], 'ro')
    # plt.xlabel('f (Hz)', fontsize=16)
    plt.ylabel('Pressure level dB SPL', fontsize=14)
    plt.title('2eme etage : h_r = 3.50m')
    plt.xlim(0, 1200)
    plt.ylim(-40, 40)
    plt.grid()

    ax = fig.add_subplot(313)
    for j in range(len(d_sr)):
        plt.plot(f[0, :], l_db_spl_subtract3[2, j, :])
    # plt.plot(f[0, :], avg_sub[:], 'ro')
    plt.xlabel('f (Hz)', fontsize=16)
    plt.ylabel('Pressure level dB SPL', fontsize=14)
    plt.title('3eme etage : h_r = 4.50m')
    plt.xlim(0, 1200)
    plt.ylim(-40, 40)
    plt.grid()

    print 10. * np.log10(np.sum(10 ** (np.abs(l_db_spl_subtract1[5, 1, :10])/10.)))

    plt.show()
