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


def levels(d_sr, d_sr_ref, h_r, case, method):
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
    p0 = np.load(os.path.join(res_path, 'p_ref.npy'))
    p1 = np.load(os.path.join(res_path, 'p_fac.npy'))

    P0, f = zip(*(basic_fft(p0[i, :], t, Ts, 0, 'hanning')
                  for i in range(len(d_sr_ref))))
    P0 = np.asarray(P0)
    l_db_spl0 = [20. * np.log10(np.sum(np.abs(P0[i, :])) / (2. * 10**-5))
                for i in range(len(d_sr_ref))]

    print p1.shape[0], p1.shape[1]

    P1 = np.zeros((p1.shape[0], p1.shape[1], p1.shape[2]-1), dtype=np.float64)
    for i in range(len(d_sr)):
        for j in range(len(h_r)):
            P1[i, j, :], f = basic_fft(p1[i, j, :], t, Ts, 0, 'hanning')
    # P1 = np.asarray(P1)
    # l_db_spl1 = [20. * np.log10(np.sum(np.abs(P1[i, :])) / (2. * 10**-5))
    #             for i in range(len(d_sr))]

    print np.shape(P0), np.shape(P1)

    f = np.asarray(f)

    frequencies = [63., 80., 100., 125., 160., 200.,
                   250., 315., 400., 500., 630., 800., 1000., 1250., 1600.,
                   2000.]
    level0_third = np.zeros((len(d_sr), len(frequencies)), dtype=np.float64)
    level1_third = np.zeros((len(d_sr), len(h_r), len(frequencies)), dtype=np.float64)
    for i in range(len(d_sr)):
        f_third, level0_third[i, :] = third_octaves(p0[i, :], 1. / Ts,
                                                    frequencies=frequencies,
                                                    ref=2e-05)
        for j in range(len(h_r)):
            f_third, level1_third[i, j, :] = third_octaves(p1[i, j, :], 1. / Ts,
                                                        frequencies=frequencies,
                                                        ref=2e-05)

    l_diff = np.zeros((len(d_sr) * len(h_r) * len(d_sr_ref) + 1, len(frequencies)),
                      dtype=np.float64)
    count = 0
    for i in range(len(d_sr)):
        for j in range(len(h_r)):
            for k in range(len(d_sr_ref)):
                l_diff[count, :] = level1_third[i, j, :] - level0_third[k, :]
                count += 1

    l_diff_sort = np.sort(l_diff, axis=0)
    l_10 = np.percentile(l_diff_sort, 10, axis=0)
    l_50 = np.percentile(l_diff_sort, 50, axis=0)
    l_90 = np.percentile(l_diff_sort, 90, axis=0)

    f_third_oct = [int(freq) for freq in frequencies]

    data = l_diff_sort
    mpl_fig = plt.figure('0')
    ax = mpl_fig.add_subplot(111)
    ax.boxplot(data)
    plt.xticks(range(1, 17), f_third_oct, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylim(-20, 20)
    plt.ylabel(r'$L_{r1-6} - L_{r0}$ (dB)', fontsize=16)
    plt.grid()
    mpl_fig.tight_layout()
    plt.savefig(os.path.join(res_path, 'l_all_rec.eps'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'l_all_rec.png'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'l_all_rec.pdf'), transparent=True, bbox_inches='tight', pad_inches=0)

    l_diff_mic_fac = np.zeros((3, len(h_r) * len(d_sr_ref) + 1, len(frequencies)),
                      dtype=np.float64)
    count = 0
    for j in range(len(h_r)):
        for k in range(len(d_sr_ref)):
            l_diff_mic_fac[0, count, :] = level1_third[0, j, :] - level0_third[k, :]
            l_diff_mic_fac[1, count, :] = level1_third[10, j, :] - level0_third[k, :]
            l_diff_mic_fac[2, count, :] = level1_third[-1, j, :] - level0_third[k, :]
            count += 1

    data1 = np.sort(l_diff_mic_fac[0, :, :], axis=0)
    data2 = np.sort(l_diff_mic_fac[1, :, :], axis=0)
    data3 = np.sort(l_diff_mic_fac[2, :, :], axis=0)
    mpl_fig = plt.figure('3')
    ax = mpl_fig.add_subplot(311)
    ax.boxplot(data1)
    plt.title(r'(a) Dist. btw. facade and receiver: $d_{fr}$ = %0.1f m.' % d_sr[0])
    plt.xticks(range(1, 17), f_third_oct, rotation='vertical')
    plt.ylim(-20, 20)
    plt.ylabel(r'$L_{r1-6} - L_{r0}$ (dB)')
    plt.grid()
    mpl_fig.tight_layout()

    ax = mpl_fig.add_subplot(312)
    ax.boxplot(data2)
    plt.title(r'(b) Dist. btw. facade and receiver: $d_{fr}$ = %0.1f m.' % d_sr[10])
    plt.xticks(range(1, 17), f_third_oct, rotation='vertical')
    plt.ylim(-20, 20)
    plt.ylabel(r'$L_{r1-6} - L_{r0}$ (dB)')
    plt.grid()
    mpl_fig.tight_layout()

    ax = mpl_fig.add_subplot(313)
    ax.boxplot(data3)
    plt.title(r'(c) Dist. btw. facade and receiver: $d_{fr}$ = %0.1f m.' % d_sr[-1])
    plt.xticks(range(1, 17), f_third_oct, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylim(-20, 20)
    plt.ylabel(r'$L_{r1-6} - L_{r0}$ (dB)')
    plt.grid()
    mpl_fig.tight_layout()
    plt.savefig(os.path.join(res_path, 'l_3dfr.eps'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'l_3dfr.png'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'l_3dfr.pdf'), transparent=True, bbox_inches='tight', pad_inches=0)


    print np.shape(level0_third), np.shape(level1_third)

    frequencies = np.asarray(frequencies)
    fig = plt.figure('1')
    ax = fig.add_subplot(311)
    plt.semilogx(frequencies-1.5,
                 np.median(level1_third[0, 1, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[0, 8, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[0, 18, :] - level0_third[:, :], axis=0))

    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    plt.legend(('h_r = %0.1f m' % h_r[1], 'h_r = %0.1f m' % h_r[8], 'h_r = %0.1f m' % h_r[18]))
    # plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.title('Distance du mic-facade = %0.2f m' % d_sr[0])
    plt.grid()

    ax = fig.add_subplot(312)
    plt.semilogx(frequencies-1.5,
                 np.median(level1_third[10, 2, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[10, 6, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[10, 20, :] - level0_third[:, :], axis=0))
    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    # plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.title('Distance du mic.-facade = %0.2f m' % d_sr[10])
    plt.grid()

    ax = fig.add_subplot(313)
    plt.semilogx(frequencies-1.5,
                 np.median(level1_third[19, 2, :]-level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[19, 6, :] - level0_third[:, :], axis=0))
    plt.semilogx(frequencies - 1.5,
                 np.median(level1_third[19, 20, :] - level0_third[:, :], axis=0))
    # ax.bar(frequencies-1.5, np.median(level1_third[0, :]-level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies,
    #        np.median(level2_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.bar(frequencies + 1.5,
    #        np.median(level3_third[0, :] - level0_third[:, :], axis=0), width=5)
    # ax.set_xscale("log")
    plt.xticks(frequencies, frequencies, rotation='vertical')
    plt.xlabel('f (Hz)', fontsize=14)
    plt.ylabel('dB (SPL)', fontsize=14)
    plt.title('Distance du mic-facade = %0.2f m' % d_sr[19])
    plt.grid()

    plt.show()
