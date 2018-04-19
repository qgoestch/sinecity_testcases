# -*- coding: utf-8 -*-
##
# \file     plot_time_signals.py
# \title    Comparison of the time signals at chosen receiver locations.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar), LAUM (Le Mans Université)
# \date     2017, 12 Oct.
##
import numpy as np
import os
import site
base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))


def plot_ts_2(p_an, d_sr, c, freq, case, h_num, num_meth):
    """

    Plot the time signals for selected receivers - used in case 1 for
    the geometrical spreading verification.
    :param p_an: analytic pressure calculated in errors_calc2_geospr.py (Pa).
    :type p_an: 2d-array of floats [distance, time]
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param c: sound speed (m.s-1).
    :type c: float
    :param freq: frequency sequence (Hz).
    :type freq: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param h_num: spatial step index.
    :type h_num: int
    :param num_meth: name of the numerical method: fdtd or tlm
    :type num_meth: string
    :return: the time signals for axial and diagonal
    """
    # =========================================================================
    # -------------------------- Check time signals ---------------------------
    # =========================================================================
    from matplotlib import pyplot as plt
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, num_meth)
    t = np.load(os.path.join(res_path, 't_%i.npy' % (h_num)))
    dist_idx = 0
    p_axial = np.load(os.path.join(res_path, 'p_axial_%i.npy' % (h_num)))
    p_diag = np.load(os.path.join(res_path, 'p_diag_%i.npy' % (h_num)))

    diff = np.zeros((p_axial.shape[0], p_axial.shape[1]))
    for j in range(p_axial.shape[0]):
        for i in range(p_axial.shape[1]):
            diff[j, i] = np.abs(p_diag[j, i] - p_axial[j, i])

    if np.abs(len(t) - len(p_axial[dist_idx, : ]))!=0:
        idx_cor = np.abs(len(t) - len(p_axial[dist_idx, : ]))
    else:
        idx_cor = 0

    dt_n = 4
    if h_num == 4 or h_num == 3 or h_num == 2:
        dt_n = 2

    fig = plt.figure('Time signals 2 %s %i ' % (num_meth, h_num))
    ax = fig.add_subplot(211)
    plt.plot(t, p_an[dist_idx, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx, idx_cor:])
    plt.legend(('Analytic', 'Diagonal'))
    plt.plot(t, p_an[dist_idx + 1, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 1, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 1, idx_cor:])
    plt.plot(t, p_an[dist_idx + 2, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 2, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 2, idx_cor:])
    plt.plot(t, p_an[dist_idx + 3, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 3, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 3, idx_cor:])
    plt.plot(t, p_an[dist_idx + 4, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 4, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 4, idx_cor:])
    plt.plot(t, p_an[dist_idx + 5, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 5, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 5, idx_cor:])
    plt.plot(t, p_an[dist_idx + 6, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 6, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 6, idx_cor:])
    plt.plot(t, p_an[dist_idx + 7, idx_cor:], 'k:')
    # plt.plot(t[::dt_n], p_axial[dist_idx + 7, :-1:dt_n], '--')
    plt.plot(t, p_diag[dist_idx + 7, idx_cor:])
    plt.ylabel('Pressure (Pa)', fontsize=12)
    plt.xlim(0.0218, 0.0335)

    ax = fig.add_subplot(212)
    plt.plot(t, p_an[dist_idx, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx, idx_cor:])
    plt.legend(('Analytic', 'Axial'))
    plt.plot(t, p_an[dist_idx + 1, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 1, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 1, idx_cor:])
    plt.plot(t, p_an[dist_idx + 2, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 2, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 2, idx_cor:])
    plt.plot(t, p_an[dist_idx + 3, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 3, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 3, idx_cor:])
    plt.plot(t, p_an[dist_idx + 4, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 4, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 4, idx_cor:])
    plt.plot(t, p_an[dist_idx + 5, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 5, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 5, idx_cor:])
    plt.plot(t, p_an[dist_idx + 6, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 6, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 6, idx_cor:])
    plt.plot(t, p_an[dist_idx + 7, idx_cor:], 'k:')
    plt.plot(t[::dt_n], p_axial[dist_idx + 7, :-1:dt_n], '--')
    # plt.plot(t, p_diag[dist_idx + 7, idx_cor:])
    plt.ylabel('Pressure (Pa)', fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    plt.xlim(0.0218, 0.0335)
    # plt.close()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'time_signals_2_%s_%i_geospr.png'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'time_signals_2_%s_%i_geospr.pdf'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'time_signals_2_%s_%i_geospr.eps'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)

    plt.show()


def plot_ts(d_sr, c, freq, case, h_num, num_meth):
    """

    Plot the time signals for selected receivers - used in case 1 for
    the geometrical spreading verification.

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param c: sound speed (m.s-1).
    :type c: float
    :param freq: frequency sequence (Hz).
    :type freq: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param h_num: spatial step index.
    :type h_num: int
    :param num_meth: name of the numerical method: fdtd or tlm
    :type num_meth: string
    :return: the time signals for axial and diagonal
    """
    # =========================================================================
    # -------------------------- Check time signals ---------------------------
    # =========================================================================
    from matplotlib import pyplot as plt
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, num_meth)
    t = np.load(os.path.join(res_path, 't_%i.npy' % (h_num)))
    Ts = np.load(os.path.join(res_path, 'Ts_%i.npy' % (h_num)))
    dist_idx = 0
    p_exact = [1. / (d_sr[dist_idx]) ** np.log2(10) *
               np.exp(-(np.pi ** 2) *
               ((freq / 2.) * t[n - int(round((d_sr[dist_idx] / c) / Ts))] - 1) ** 2)
               for n in range(len(t))]
    p_axial = np.load(os.path.join(res_path, 'p_axial_%i.npy' % (h_num)))
    p_diag = np.load(os.path.join(res_path, 'p_diag_%i.npy' % (h_num)))

    diff = np.zeros((p_axial.shape[0], p_axial.shape[1]))
    for j in range(p_axial.shape[0]):
        for i in range(p_axial.shape[1]):
            diff[j, i] = np.abs(p_diag[j, i] - p_axial[j, i])

    if np.abs(len(t) - len(p_axial[dist_idx, : ]))!=0:
        idx_cor = np.abs(len(t) - len(p_axial[dist_idx, : ]))
    else:
        idx_cor = 0

    dt_n = 4
    if h_num == 4 or h_num == 3 or h_num == 2:
        dt_n = 2

    fig = plt.figure('Time signals %s %i ' % (num_meth, h_num))
    ax = fig.add_subplot(211)
    plt.plot(t[::dt_n], p_axial[dist_idx, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx, idx_cor:])
    plt.legend(('Axial', 'Diagonal'))
    plt.plot(t[::dt_n], p_axial[dist_idx + 1, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 1, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 2, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 2, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 3, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 3, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 4, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 4, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 5, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 5, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 6, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 6, idx_cor:])
    plt.plot(t[::dt_n], p_axial[dist_idx + 7, :-1:dt_n], 'k*')
    plt.plot(t, p_diag[dist_idx + 7, idx_cor:])
    plt.ylabel('Pressure (Pa)', fontsize=12)

    ax = fig.add_subplot(212)
    plt.plot(t, diff[dist_idx, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 1, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 2, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 3, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 4, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 5, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 6, idx_cor:], '--')
    plt.plot(t, diff[dist_idx + 7, idx_cor:], '--')
    plt.ylim(0, 0.1)
    plt.ylabel(r'|Difference| (Pa)', fontsize=12)
    plt.xlabel('Time (s)', fontsize=12)
    # plt.close()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'time_signals_%s_%i_geospr.png'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'time_signals_%s_%i_geospr.pdf'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'time_signals_%s_%i_geospr.eps'
                             % (num_meth, h_num)), transparent=True,
                bbox_inches='tight', pad_inches=0)

    plt.show()


def plot_ts_basic(t, p_1, p_2):
    """

    :param t:
    :type t:
    :param p_1:
    :type p_1:
    :param p_2:
    :type p_2:
    :param save_fig:
    :type save_fig:
    :return:
    :rtype:
    """
    from matplotlib import pyplot as plt
    fig = plt.figure('Time signals')
    plt.plot(t, p_1)
    plt.plot(t, p_2)
    plt.ylabel('Pressure (Pa)', fontsize=12)


    plt.show()