# -*- coding: utf-8 -*-
##
# \file     plot_att_spect.py
# \title    Save the attenuation spectrums at 6 positions
#           for a given grid in case3: ground reflection.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 17 Oct.
##

import numpy as np
from matplotlib import pyplot as plt
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def attenuation_spectrums(f_fdtd, P_T_fdtd, P_F_fdtd, f_tlm, P_T_tlm, P_F_tlm,
                          f_an, P_T_an, P_F_an, case):
    """
    Plot the attenuation spectrums for case 3 in order to compare the numerical
    methods to the analytic solution. Six receiver positions are selected,
    which results in 6 subplots where the spectrums might be overlapped.

    :param f_fdtd: frequency sequence for the fdtd (Hz).
    :type f_fdtd: 1d list of floats
    :param P_T_fdtd: total pressure in presence of BC for the fdtd (Pa).
    :type P_T_fdtd: 1d list [n] of floats
    :param P_F_fdtd: free-field pressure for the fdtd (Pa).
    :type P_F_fdtd: 1d list [n] of floats
    :param f_tlm: frequency sequence for the tlm (Hz).
    :type f_tlm: 1d list of floats
    :param P_T_tlm: total pressure in presence of BC for the tlm (Pa).
    :type P_T_tlm: 1d list [n] of floats
    :param P_F_tlm: free-field pressure for the tlm (Pa).
    :type P_F_tlm: 1d list [n] of floats
    :param f_an: frequency sequence for the analytic solution (Hz).
    :type f_an: 1d list of floats
    :param P_T_an: total field in presence of BC for the analytic sol. (Pa).
    :type P_T_an: 1d list [n] of floats
    :param P_F_an: free-field pressure for the analytic solution (Pa).
    :type P_F_an: 1d list [n] of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :return: the attenuation spectrums relative to free-field propagation.
    """
    # =========================================================================
    #   Plot the Attenuation spectrums (rel. to free field) for all methods
    # =========================================================================
    ATT_fdtd = 10. * np.log10(np.abs(P_T_fdtd / P_F_fdtd)**2)
    ATT_tlm = 10. * np.log10(np.abs(P_T_tlm / P_F_tlm)**2)
    ATT_an = 10. * np.log10(np.abs(P_T_an / P_F_an)**2)

    print np.shape(ATT_fdtd)
    d_idx_1 = 4
    h_idx_1 = 3
    d_idx_2 = 9
    h_idx_2 = h_idx_1
    d_idx_3 = 15
    h_idx_3 = h_idx_1
    d_idx_4 = 4
    h_idx_4 = 4
    d_idx_5 = 9
    h_idx_5 = h_idx_4
    d_idx_6 = 15
    h_idx_6 = h_idx_4
    fig = plt.figure('Attenuation spectrums', figsize=(14, 6))
    ax = fig.add_subplot(231)
    ax.plot(f_an, ATT_an[d_idx_1, h_idx_1, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_1, h_idx_1, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_1, h_idx_1, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.legend(('Analytic', 'TLM', 'FDTD'), fontsize=14)
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(232)
    ax.plot(f_an, ATT_an[d_idx_2, h_idx_2, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_2, h_idx_2, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_2, h_idx_2, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(233)
    ax.plot(f_an, ATT_an[d_idx_3, h_idx_3, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_3, h_idx_3, :], 'rs',markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_3, h_idx_3, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(234)
    ax.plot(f_an, ATT_an[d_idx_4, h_idx_4, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_4, h_idx_4, :], 'rs',markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_4, h_idx_4, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(235)
    ax.plot(f_an, ATT_an[d_idx_5, h_idx_5, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_5, h_idx_5, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_5, h_idx_5, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    ax = fig.add_subplot(236)
    ax.plot(f_an, ATT_an[d_idx_6, h_idx_6, :])
    ax.plot(f_tlm, ATT_tlm[d_idx_6, h_idx_6, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(f_fdtd, ATT_fdtd[d_idx_6, h_idx_6, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Attenuation (dB)", fontsize=12)
    plt.grid()
    plt.ylim(-40, 10)
    plt.tight_layout()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'att_spect_ground.eps'),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'att_spect_ground.png'),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'att_spect_ground.pdf'),
                transparent=True, bbox_inches='tight', pad_inches=0)

    plt.show()
