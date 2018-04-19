# -*- coding: utf-8 -*-
##
# \file     errors_calc_ground.py
# \title    Calculation of the errors and norms for the case3: ground reflection.
# \author   Pierre Chobeau
# \version  0.2
# \date     2017, 09 Aug.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                             'num_methods', 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_ground_arrays

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from fft_td_sig_arrays import fft_conv_arrays
from miki_imp_model import Miki
from error_norm_freq import error, two_norm, max_norm
from obs_ord_acc_freq import obs_ord_acc

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_time_signals import plot_ts
from plot_errors_norms_cfa2018 import plot_errors_norms_cfa2018
from plot_att_spect import attenuation_spectrums
from plot_obs_ord_acc_freq import plot_obs_ord_acc_freq


def error_calc(d_sr, h_s, h_r, h_set, rho, c, sigma, case,
               disp_att_spect, disp_errors):
    """
    Calculation of the error (avaraged, two-norm and max-norm) for case 3.
    The error is given as a function of frequency for each grid and each
    method.
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param sigma: pecific airflow resistivity (kNm-4s==CGS).
    :type sigma: float
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param disp_att_spect: display the attenuation spectrums at 6 selected locations.
    :type disp_att_spect: bool
    :param disp_errors: display the errors in norms.
    :type disp_errors: bool

    :param P_T_fdtd: total pressure in presence of BC for the fdtd,
    (3+1)d list [h_set][d_sr,h_r,n].
    :param P_F_fdtd: free-field pressure for the fdtd,
    (3+1)d list [h_set][d_sr,h_r,n].
    :param f_fdtd: frequency sequence for the fdtd, 1d list.
    :param P_T_tlm: total pressure in presence of BC for the tlm,
    (3+1)d list [h_set][d_sr,h_r,n].
    :param P_F_tlm: free-field pressure for the tlm,
    (3+1)d list [h_set][d_sr,h_r,n].
    :param f_tlm: frequency sequence for the tlm, 1d list.

    :param max_freq_idx: index of the maximal frequency in relation with the
    source and the grid max size, scalar.

    :param omega: angular frequency, 1d-list.
    :param k_f: wave number, 1d-list.
    :param Zg: surface impedance from Miki's model, 1d-list.
    :param k_m: wave number from Miki's model, 1d list.
    :param P_T_an: total pressure in presence of BC - analytic solution,
    (3+1)d list [h_set][d_sr,h_r,n].
    :param P_F_an: free-field pressure - analytic solution,
    (3+1)d list [h_set][d_sr,h_r,n].

    :param avg_error_tlm: error averaged over all receivers for the TLM for
    each spatial step, 1d-array.
    :param avg_error_fdtd: error averaged over all receivers for the FDTD for
    each spatial step, 1d-array.

    :param two_norm_tlm: relative error in the 2-norm for the TLM as function
                            of frequency for each spatial step, 1d-array.
    :param two_norm_fdtd: relative error in the 2-norm for the FDTD as function
                            of frequency for each spatial step, 1d-array.
    :param max_norm_tlm: relative error in the max-norm for the TLM as function
                            of frequency for each spatial step, 1d-array.
    :param max_norm_fdtd: relative error in the max-norm for the FDTD as function
                            of frequency for each spatial step, 1d-array.

    :param ord_acc_tlm: order of accuracy between two consecutive grids in
    2-norm for the TLM as function of frequency for each spatial step,
    1d-array.
    :param ord_acc_fdtd: order of accuracy between two consecutive grids in
    2-norm for the FDTD as function of frequency for each spatial step,
    1d-array.
    :param ord_acc_tlm_max: order of accuracy between two consecutive grids in
    max-norm for the TLM as function of frequency for each spatial step,
    1d-array.
    :param ord_acc_fdtd_max: order of accuracy between two consecutive grids in
    max-norm for the FDTD as function of frequency for each spatial step,
    1d-array.
    :param freq_idx: chosen frequency index within the range(max_freq_idx),
    float.

    :return The log-log plot of the two-norm and max-norm
    (via plot_errors_norms.py), and the respective orders of accuracy for the
    chosen frequency. It shows also the attenuation spectrums for 6 receiver
    locations (via plot_att_spect.py).
    """
    #   Load the numerical results
    P_T_fdtd, P_F_fdtd, f_fdtd = zip(*(fft_conv_arrays('fdtd', case, i)
                                       for i in range(len(h_set))))
    P_T_tlm, P_F_tlm, f_tlm = zip(*(fft_conv_arrays('tlm', case, i)
                                    for i in range(len(h_set))))

    #   Setting the frequency span used for all grids
    max_freq_idx = 40
    f_fdtd = f_fdtd[0][:max_freq_idx]
    f_tlm = f_tlm[0][:max_freq_idx]
    print f_tlm
    #   Calculation of the analytic solutions
    omega = 2. * np.pi * f_tlm
    k_f = omega / c
    Zg, k_m = Miki(-1, f_tlm, sigma, rho, c)
    P_F_an, P_T_an = analytic_solution_ground_arrays(d_sr, h_s, h_r,
                                                     Zg / (rho * c), k_f)

    #   Errors and norms
    one_norm_tlm = np.zeros((len(range(max_freq_idx)), len(h_set)))
    one_norm_fdtd = np.zeros((len(range(max_freq_idx)), len(h_set)))
    two_norm_tlm = np.zeros((len(range(max_freq_idx)), len(h_set)))
    two_norm_fdtd = np.zeros((len(range(max_freq_idx)), len(h_set)))
    max_norm_tlm = np.zeros((len(range(max_freq_idx)), len(h_set)))
    max_norm_fdtd = np.zeros((len(range(max_freq_idx)), len(h_set)))
    ord_acc_tlm_one = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    ord_acc_fdtd_one = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    ord_acc_tlm_two = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    ord_acc_fdtd_two = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    ord_acc_tlm_max = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    ord_acc_fdtd_max = np.zeros((len(range(max_freq_idx)), len(h_set) - 1))
    # =========================================================================
    #   Calculation of the errors and norms using numpy.linalg.norm
    # Work for the 2-norm, not for the max-norm that is only doable for
    # vectors (not possible for the arrays >= 2D)
    # =========================================================================
    for l in range(len(h_set)):
        for i in range(max_freq_idx):
            one_norm_fdtd[i, l] = np.linalg.norm((P_T_fdtd[l][:, :, i] -
                                                      P_T_an[:, :, i]) *
                                                      h_set[l] ** 2, ord=1)
            one_norm_tlm[i, l] = np.linalg.norm((P_T_tlm[l][:, :, i] -
                                                      P_T_an[:, :, i]) *
                                                      h_set[l] ** 2, ord=1)

            two_norm_fdtd[i, l] = np.linalg.norm((P_T_fdtd[l][:, :, i] -
                                                      P_T_an[:, :, i]) *
                                                      h_set[l] ** 2, ord=2)
            two_norm_tlm[i, l] = np.linalg.norm((P_T_tlm[l][:, :, i] -
                                                      P_T_an[:, :, i]) *
                                                      h_set[l] ** 2, ord=2)

            max_norm_fdtd[i, l] = np.linalg.norm((P_T_fdtd[l][:, :, i] -
                                                  P_T_an[:, :, i]) *
                                                 h_set[l] ** 2, ord=np.inf)
            max_norm_tlm[i, l] = np.linalg.norm((P_T_tlm[l][:, :, i] -
                                                 P_T_an[:, :, i]) *
                                                h_set[l] ** 2, ord=np.inf)
    # =========================================================================
    #   Calculation of the order of accuracy
    # =========================================================================
    for i in range(max_freq_idx):
        for l in range(len(h_set)-1):
            ord_acc_tlm_one[i, l] = np.log(one_norm_tlm[i, l + 1] /
                                           one_norm_tlm[i, l]) / \
                                    np.log(h_set[l + 1]/h_set[l])
            ord_acc_fdtd_one[i, l] = np.log(one_norm_fdtd[i, l + 1] /
                                            one_norm_fdtd[i, l]) / \
                                     np.log(h_set[l + 1]/h_set[l])
            ord_acc_tlm_two[i, l] = np.log(two_norm_tlm[i, l + 1] /
                                           two_norm_tlm[i, l]) / \
                                    np.log(h_set[l + 1]/h_set[l])
            ord_acc_fdtd_two[i, l] = np.log(two_norm_fdtd[i, l + 1] /
                                            two_norm_fdtd[i, l]) / \
                                     np.log(h_set[l + 1]/h_set[l])
            ord_acc_tlm_max[i, l] = np.log(max_norm_tlm[i, l + 1] /
                                           max_norm_tlm[i, l]) / \
                                    np.log(h_set[l + 1]/h_set[l])
            ord_acc_fdtd_max[i, l] = np.log(max_norm_fdtd[i, l + 1] /
                                            max_norm_fdtd[i, l]) / \
                                     np.log(h_set[l + 1]/h_set[l])
    # =========================================================================
    #   Observed order of accuracy as a function of frequency
    # =========================================================================
    p_obs_fdtd = obs_ord_acc(P_T_fdtd, f_fdtd)
    p_obs_tlm = obs_ord_acc(P_T_tlm, f_tlm)
    print np.shape(p_obs_fdtd)
    plot_obs_ord_acc_freq(p_obs_fdtd, p_obs_tlm, f_fdtd, case)

    freq_idx = 4
    print 'The frequency considered for the error plots is f = %.2f Hz' \
          % f_fdtd[freq_idx]

    if disp_errors:
        plot_errors_norms_cfa2018(h_set, np.average(one_norm_tlm, axis=0),
                          np.average(one_norm_fdtd, axis=0),
                          np.average(two_norm_tlm, axis=0),
                          np.average(two_norm_fdtd, axis=0),
                          np.average(max_norm_tlm, axis=0),
                          np.average(max_norm_fdtd, axis=0),
                          np.average(ord_acc_tlm_one, axis=0),
                          np.average(ord_acc_fdtd_one, axis=0),
                          np.average(ord_acc_tlm_two, axis=0),
                          np.average(ord_acc_fdtd_two, axis=0),
                          np.average(ord_acc_tlm_max, axis=0),
                          np.average(ord_acc_fdtd_max, axis=0), case)

    if disp_att_spect:
        grid_idx = 0
        attenuation_spectrums(f_fdtd, P_T_fdtd[grid_idx][:, :, :max_freq_idx],
                              P_F_fdtd[grid_idx][:, :, :max_freq_idx],
                              f_tlm, P_T_tlm[grid_idx][:, :, :max_freq_idx],
                              P_F_tlm[grid_idx][:, :, :max_freq_idx],
                              f_tlm, P_T_an[:, :, :], P_F_an[:, :, :], case)
