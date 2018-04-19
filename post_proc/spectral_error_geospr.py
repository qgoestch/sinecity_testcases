# -*- coding: utf-8 -*-
##
# \file     spectral_error_geospr.py
# \title    Calculation of the errors as a function of frequency for
#           case1: geometrical spreading.
# \author   Pierre Chobeau
# \version  0.2
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2018, 12 Jan.
##
import numpy as np
import os
import site
from matplotlib import pyplot as plt

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                             'num_methods', 'analytic')

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from fft_td_sig_arrays import basic_fft, amplitude_spectrum

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from plot_time_signals import plot_ts


def spectral_error(d_sr, case, f_max_src):
    """

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param f_max_src: approximated maximale frequency of the source signal (Hz).
    :type f_max_src: float
    :return:
    :rtype:
    """

    # =========================================================================
    #   Load the time signals for each method
    # =========================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, 'fdtd')
    l = 0
    p_axial_fdtd = np.load(os.path.join(res_path, 'p_axial_%i.npy' % l))
    p_diag_fdtd = np.load(os.path.join(res_path, 'p_diag_%i.npy' % l))
    t = np.load(os.path.join(res_path, 't_%i.npy' % l))
    Ts = np.load(os.path.join(res_path, 'Ts_%i.npy' % l))

    print len(t), np.shape(p_axial_fdtd)

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, 'tlm')
    p_axial_tlm = np.load(os.path.join(res_path, 'p_axial_%i.npy' % l))
    p_diag_tlm = np.load(os.path.join(res_path, 'p_diag_%i.npy' % l))
    t = np.load(os.path.join(res_path, 't_%i.npy' % l))
    Ts = np.load(os.path.join(res_path, 'Ts_%i.npy' % l))

    print len(t), np.shape(p_axial_tlm)

    # =========================================================================
    #   Calculation of the FFTs
    # =========================================================================
    f_test, P_test = amplitude_spectrum(abs(p_axial_fdtd[0, :]), 1/Ts,
                                        p_diag_fdtd.shape[1])
    P_test = np.abs(P_test/len(p_axial_fdtd[0, :]))

    P_axial_fdtd, f_ax_fdtd = zip(*(basic_fft(p_axial_fdtd[i, :], t, Ts, 0, 'hanning')
                            for i in range(len(d_sr))))
    P_diag_fdtd, f_di_fdtd = zip(*(basic_fft(p_diag_fdtd[i, :], t, Ts, 0, 'hanning')
                            for i in range(len(d_sr))))
    P_axial_tlm, f_ax_tlm = zip(*(basic_fft(p_axial_tlm[i, :], t, Ts, 0, 'hanning')
                            for i in range(len(d_sr))))
    P_diag_tlm, f_di_tlm = zip(*(basic_fft(p_diag_tlm[i, :], t, Ts, 0, 'hanning')
                            for i in range(len(d_sr))))

    f_di_tlm = np.asarray(f_di_tlm)
    P_axial_fdtd = np.asarray(P_axial_fdtd)
    print np.shape(P_axial_fdtd), np.shape(f_di_tlm)
    P_diag_fdtd = np.asarray(P_diag_fdtd)
    P_axial_tlm = np.asarray(P_axial_tlm)
    P_diag_tlm = np.asarray(P_diag_tlm)

    #   Setting the frequency span used for all grids
    max_freq_idx = int(f_di_tlm.shape[1]/2.)
    print max_freq_idx
    f = f_di_tlm[0, :max_freq_idx]
    # print f

    LdB_rel_axial_fdtd = np.zeros((max_freq_idx))
    LdB_rel_axial_tlm = np.zeros((max_freq_idx))
    LdB_rel_diag_fdtd = np.zeros((max_freq_idx))
    LdB_rel_diag_tlm = np.zeros((max_freq_idx))
    LdB_rel_axial_fdtd_abs = np.zeros((P_axial_fdtd.shape[0], max_freq_idx))
    LdB_rel_axial_tlm_abs = np.zeros((P_axial_fdtd.shape[0], max_freq_idx))
    LdB_rel_diag_fdtd_abs = np.zeros((P_axial_fdtd.shape[0], max_freq_idx))
    LdB_rel_diag_tlm_abs = np.zeros((P_axial_fdtd.shape[0], max_freq_idx))
    for i in range(max_freq_idx):
        LdB_rel_axial_fdtd[i] = 10. * np.log10(P_axial_fdtd[7, i]/ \
                                                P_axial_fdtd[0, i]) ** 2
        LdB_rel_axial_tlm[i] = 10. * np.log10(P_axial_tlm[7, i]/ \
                                                P_axial_tlm[0, i]) ** 2
        LdB_rel_diag_fdtd[i] = 10. * np.log10(P_diag_fdtd[7, i]/ \
                                                P_diag_fdtd[0, i]) ** 2
        LdB_rel_diag_tlm[i] = 10. * np.log10(P_diag_tlm[7, i]/ \
                                                P_diag_tlm[0, i]) ** 2
        for j in range(P_axial_fdtd.shape[0]):
            LdB_rel_axial_fdtd_abs[j, i] = 10. * np.log10(P_axial_fdtd[j, i]/ \
                                                       (2. * 10**-5)) ** 2
            LdB_rel_axial_tlm_abs[j, i] = 10. * np.log10(P_axial_tlm[j, i]/ \
                                                      (2. * 10 ** -5)) ** 2
            LdB_rel_diag_fdtd_abs[j, i] = 10. * np.log10(P_diag_fdtd[j, i]/ \
                                                      (2. * 10 ** -5)) ** 2
            LdB_rel_diag_tlm_abs[j, i] = 10. * np.log10(P_diag_tlm[j, i]/ \
                                                     (2. * 10 ** -5)) ** 2

    print np.shape(LdB_rel_axial_fdtd)



    fig = plt.figure('Spectrums rel.', figsize=(8, 6))
    ax = fig.add_subplot(211)
    ax.plot(f, LdB_rel_axial_fdtd)
    ax.plot(f, LdB_rel_axial_tlm)
    # ax.plot(f, P_axial_tlm[0, :max_freq_idx], '--', markersize=10,
    #         markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both",ls=":")
    plt.legend(('FDTD', 'TLM'), fontsize=14)
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.xlim(0, 1. * f_max_src)
    # plt.xlim(0, 0.25*1./Ts)
    # plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(212)
    ax.plot(f, LdB_rel_diag_fdtd)
    ax.plot(f, LdB_rel_diag_tlm)
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.xlim(0, 1. * f_max_src)
    # plt.xlim(0, 0.25 * 1. / Ts)
    # plt.ylim(-40, 10)
    plt.tight_layout()

    fig = plt.figure('Spectrums abs.', figsize=(8, 6))
    ax = fig.add_subplot(211)
    ax.plot(f_test, P_test)
    # for j in range(P_axial_fdtd.shape[0] - 6):
    #     # ax.plot(f, LdB_rel_axial_fdtd_abs[j, :])
    #     ax.plot(f, P_axial_fdtd[j, :max_freq_idx])
    ax.grid(True, which="both",ls=":")
    # plt.legend(('FDTD', 'TLM'), fontsize=14)
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    # plt.xlim(0, 1. * f_max_src)
    plt.tight_layout()

    ax = fig.add_subplot(212)
    for j in range(P_axial_fdtd.shape[0]):
        ax.plot(f, P_diag_fdtd[j, :max_freq_idx])
        # ax.plot(f, LdB_rel_diag_tlm_abs[j, :])
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    # plt.xlim(0, 1. * f_max_src)
    plt.tight_layout()

    plot_ts(d_sr, 340, f_max_src, case, 0, 'fdtd')
