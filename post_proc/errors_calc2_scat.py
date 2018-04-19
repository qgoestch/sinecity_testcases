# -*- coding: utf-8 -*-
##
# \file     errors_calc2_scat.py
# \title    Calculation of the errors and norms part 2.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 24 Jan.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_scattered_pressure
from fft_td_sig_arrays import fft_conv_arrays, basic_fft
from error_norm_freq import error, two_norm, max_norm

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from plot_polar_scat import polar_plots
from plot_errors_norms_cfa2018 import plot_errors_norms_cfa2018


def error_calc2(h_set, rho, c, radius, case):
    """
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param radius: radius of the scatterer (m)
    :type radius: float
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :return:
    :rtype:
    """
    f_idx_start = 31
    f_idx_end = 40  # 750Hz ; 120 # 1500 Hz
    f_idx_end_err = f_idx_start + 2  # !!! Warning !!! Specific frequency range
                                     # for the error calculation
    one_norm_tlm = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                              len(h_set)), dtype=np.float64)
    one_norm_fdtd = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                               len(h_set)), dtype=np.float64)
    two_norm_fdtd = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                              len(h_set)), dtype=np.float64)
    max_norm_fdtd = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                              len(h_set)), dtype=np.float64)
    two_norm_tlm = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                             len(h_set)), dtype=np.float64)
    max_norm_tlm = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                             len(h_set)), dtype=np.float64)
    ord_acc_tlm_one = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                len(h_set)-1), dtype=np.float64)
    ord_acc_fdtd_one = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                 len(h_set)-1), dtype=np.float64)
    ord_acc_tlm_two = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                len(h_set)-1), dtype=np.float64)
    ord_acc_fdtd_max = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                 len(h_set)-1), dtype=np.float64)
    ord_acc_fdtd_two = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                 len(h_set)-1), dtype=np.float64)
    ord_acc_tlm_max = np.zeros((int(f_idx_end_err - f_idx_start + 1),
                                len(h_set)-1), dtype=np.float64)

    for l, h_val in enumerate(h_set[:]):
        print 'h = %f ; %i/%i' % (h_val, l,len(h_set)-1)
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'analytic', 'p_h%i' % l)
        pan = np.load(os.path.join(res_path, 'p.npy'))

        res_path_fdtd = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                     'case%i' % case, 'fdtd', 'p_h%i' % l)
        dist_rcp_grid = np.load(
            os.path.join(res_path_fdtd, 'rcvdist.npy'))
        phi_rcp_grid = np.load(
            os.path.join(res_path_fdtd, 'rcvphi.npy'))
        P_fdtd = np.load(os.path.join(res_path_fdtd, 'P_fdtd.npy'))
        f_fdtd = np.load(os.path.join(res_path_fdtd, 'f_fdtd.npy'))

        res_path_tlm = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                    'case%i' % case, 'tlm', 'p_h%i' % l)
        P_tlm = np.load(os.path.join(res_path_tlm, 'P_tlm.npy'))
        f_tlm = np.load(os.path.join(res_path_tlm, 'f_tlm.npy'))

        norma_angle_idx = int(pan.shape[0] / 2)
        norma_dist_idx =  int(pan.shape[1] / 2)
        norma_freq_idx = int((f_idx_end - f_idx_start)/2. + f_idx_start)

        print f_fdtd[f_idx_start: f_idx_end_err]
        print f_tlm[norma_freq_idx], f_fdtd[norma_freq_idx]
        # print 'distance idx = %i' % norma_dist_idx

        Mag_Norm_fdtd = np.abs(pan[norma_angle_idx, norma_dist_idx, norma_freq_idx]) / np.abs(
            P_fdtd[norma_angle_idx, norma_dist_idx, norma_freq_idx])
        Mag_Norm_tlm = np.abs(pan[norma_angle_idx, norma_dist_idx, norma_freq_idx]) / np.abs(
            P_tlm[norma_angle_idx, norma_dist_idx, norma_freq_idx])

        polar_plots_on = False # if true: set f_idx_start = 1 & f_idx_end = 40
        if l == 0:
            if polar_plots_on:
                polar_plots(c, radius, f_tlm, phi_rcp_grid, P_tlm, P_fdtd, pan,
                            f_idx_start, f_idx_end, Mag_Norm_tlm, Mag_Norm_fdtd,
                            case, 'TLM_FDTD')

        # Calculation of the errors and norms
        iangle = int(P_fdtd.shape[0] / 2)

        error_fdtd = np.zeros((iangle, len(range(norma_dist_idx, norma_dist_idx + 1)),
                               int(f_idx_end_err - f_idx_start + 1)), dtype=np.float64)
        error_tlm = np.zeros((iangle, len(range(norma_dist_idx, norma_dist_idx + 1)),
                              int(f_idx_end_err - f_idx_start + 1)), dtype=np.float64)
        for i, f_value in enumerate(f_fdtd[f_idx_start:f_idx_end_err]):
            # print f_tlm[i + f_idx_start], f_fdtd[i + f_idx_start]
            for j_idx, j in enumerate(range(norma_dist_idx, norma_dist_idx + 1)):
                for k in range(iangle):
                    error_fdtd[k, j_idx, i] = np.abs(
                        Mag_Norm_fdtd * np.abs(P_fdtd[k, j, i]) -
                        np.abs(pan[k, j, i])) ** 2
                    error_tlm[k, j_idx, i] = np.abs(
                        Mag_Norm_tlm * np.abs(P_tlm[k, j, i]) -
                        np.abs(pan[k, j, i])) ** 2
            one_norm_fdtd[i, l] = np.linalg.norm((error_fdtd[:, :, i]) *
                                                      h_set[l] ** 2, ord=1)
            one_norm_tlm[i, l] = np.linalg.norm((error_tlm[:, :, i]) *
                                                      h_set[l] ** 2, ord=1)

            two_norm_fdtd[i, l] = np.linalg.norm((error_fdtd[:, :, i]) *
                                                      h_set[l] ** 2, ord=2)
            two_norm_tlm[i, l] = np.linalg.norm((error_tlm[:, :, i]) *
                                                      h_set[l] ** 2, ord=2)

            max_norm_fdtd[i, l] = np.linalg.norm((error_fdtd[:, :, i]) *
                                                 h_set[l] ** 2, ord=np.inf)
            max_norm_tlm[i, l] = np.linalg.norm((error_tlm[:, :, i]) *
                                                h_set[l] ** 2, ord=np.inf)

    for i, f_value in enumerate(f_fdtd[f_idx_start:f_idx_end_err]):
        for l in range(len(h_set) - 1):
            ord_acc_tlm_one[i, l] = np.log(one_norm_tlm[i, l + 1] /
                                           one_norm_tlm[i, l]) / \
                                    np.log(h_set[l + 1]/h_set[l])
            ord_acc_fdtd_one[i, l] = np.log(one_norm_fdtd[i, l + 1] /
                                            one_norm_fdtd[i, l]) / \
                                     np.log(h_set[l + 1]/h_set[l])

            ord_acc_tlm_two[i, l] = np.log(
                two_norm_tlm[i, l + 1] / two_norm_tlm[i, l]) / np.log(
                h_set[l + 1] / h_set[l])
            ord_acc_fdtd_two[i, l] = np.log(
                two_norm_fdtd[i, l + 1] / two_norm_fdtd[i, l]) / np.log(
                h_set[l + 1] / h_set[l])
            ord_acc_tlm_max[i, l] = np.log(
                max_norm_tlm[i, l + 1] / max_norm_tlm[i, l]) / np.log(
                h_set[l + 1] / h_set[l])
            ord_acc_fdtd_max[i, l] = np.log(
                max_norm_fdtd[i, l + 1] / max_norm_fdtd[i, l]) / np.log(
                h_set[l + 1] / h_set[l])

    print np.mean(ord_acc_tlm_two), np.mean(ord_acc_tlm_max)
    print np.mean(ord_acc_fdtd_two), np.mean(ord_acc_fdtd_max)

    freq_idx = 10
    print 'The frequency considered for the error plots is f = %.2f Hz' \
          % f_fdtd[norma_freq_idx]
    # plot_errors_norms(h_set, avg_error_tlm[norma_freq_idx, :],
    #                   avg_error_fdtd[norma_freq_idx, :],
    #                   two_norm_tlm[norma_freq_idx, :], two_norm_fdtd[norma_freq_idx, :],
    #                   max_norm_tlm[norma_freq_idx, :], max_norm_fdtd[norma_freq_idx, :],
    #                   ord_acc_tlm_two[norma_freq_idx, :], ord_acc_fdtd_two[norma_freq_idx, :],
    #                   ord_acc_tlm_max[norma_freq_idx, :], ord_acc_fdtd_max[norma_freq_idx, :],
    #                   case)

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

    polar_plots(c, radius, f_tlm, phi_rcp_grid, P_tlm, P_fdtd, pan, f_idx_start,
                f_idx_end, Mag_Norm_tlm, Mag_Norm_fdtd, case, 'TLM_FDTD')
