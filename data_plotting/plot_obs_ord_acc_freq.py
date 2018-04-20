# -*- coding: utf-8 -*-
##
# \file     plot_obs_ord_acc_freq.py
# \title    Show the observed order of accuracy as a function of frequency.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2018, 21 Mar.
##

import numpy as np
from matplotlib import pyplot as plt
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def plot_obs_ord_acc_freq_geospr(p_obs_fdtd, p_obs_tlm, freq, case):
    """
    Show the observed order of accuracy as a function of frequency.
    :param p_obs_fdtd: observed order of accuracy from the FDTD method (Pa).
    :type p_obs_fdtd: 2d-array of floats
    :param p_obs_tlm: observed order of accuracy from the TLM method (Pa).
    :type p_obs_tlm: 2d-array of floats
    :param freq: frequency sequence from the FFT (Hz).
    :type freq: list of floats
    """
    fig = plt.figure('Observed order of accuracy', figsize=(14, 6))
    ax = fig.add_subplot(231)
    ax.plot(freq, p_obs_tlm[0, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[0, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.legend(('TLM', 'FDTD'), fontsize=14)
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(232)
    ax.plot(freq, p_obs_tlm[2, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[2, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(233)
    ax.plot(freq, p_obs_tlm[4, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[4, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(234)
    ax.plot(freq, p_obs_tlm[6, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[6, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(235)
    ax.plot(freq, p_obs_tlm[8, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[8, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(235)
    ax.plot(freq, p_obs_tlm[8, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[8, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()

    ax = fig.add_subplot(236)
    ax.plot(freq, p_obs_tlm[10, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[10, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    # plt.ylim(0, 6)
    plt.tight_layout()
    plt.show()


def plot_obs_ord_acc_freq(p_obs_fdtd, p_obs_tlm, freq, case):
    """
    Show the observed order of accuracy as a function of frequency.
    :param p_obs_fdtd: observed order of accuracy from the FDTD method (Pa).
    :type p_obs_fdtd: 2d-array of floats
    :param p_obs_tlm: observed order of accuracy from the TLM method (Pa).
    :type p_obs_tlm: 2d-array of floats
    :param freq: frequency sequence from the FFT (Hz).
    :type freq: list of floats
    """

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

    fig = plt.figure('Observed order of accuracy', figsize=(14, 6))
    ax = fig.add_subplot(231)
    ax.plot(freq, p_obs_tlm[d_idx_1, h_idx_1, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_1, h_idx_1, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True,which="both",ls=":")
    plt.legend(('TLM', 'FDTD'), fontsize=14)
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    ax = fig.add_subplot(232)
    ax.plot(freq, p_obs_tlm[d_idx_2, h_idx_2, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_2, h_idx_2, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    ax = fig.add_subplot(233)
    ax.plot(freq, p_obs_tlm[d_idx_3, h_idx_3, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_3, h_idx_3, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    ax = fig.add_subplot(234)
    ax.plot(freq, p_obs_tlm[d_idx_4, h_idx_4, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_4, h_idx_4, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    ax = fig.add_subplot(235)
    ax.plot(freq, p_obs_tlm[d_idx_5, h_idx_5, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_5, h_idx_5, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    ax = fig.add_subplot(236)
    ax.plot(freq, p_obs_tlm[d_idx_6, h_idx_6, :], 'rs', markersize=5,
            markeredgewidth=1.5, markeredgecolor='r', markerfacecolor='None')
    ax.plot(freq, p_obs_fdtd[d_idx_6, h_idx_6, :], 'k+', markersize=10,
            markeredgewidth=1.5, markeredgecolor='k', markerfacecolor='None')
    ax.grid(True, which="both", ls=":")
    plt.grid()
    plt.xlabel("Frequency (Hz)", fontsize=12)
    plt.ylabel("Obs. ord. acc.", fontsize=12)
    plt.grid()
    plt.ylim(0, 4)
    plt.tight_layout()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'obs_ord_acc.eps'),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'obs_ord_acc.png'),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'obs_ord_acc.pdf'),
                transparent=True, bbox_inches='tight', pad_inches=0)

    plt.show()
