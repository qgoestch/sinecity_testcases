# -*- coding: utf-8 -*-
##
# \file     errors_calc_scat.py
# \title    Calculation of the errors and norms for the case4: scattering.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 11 Sep.
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
from plot_errors_norms import plot_errors_norms


def error_calc(h_set, rho, c, radius, case):
    """
    
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param radius: radius of the scatterer (m)
    :type radius: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    
    :param  res_path    path where to fetch the results for each method, string chain.
    :param  p_t_fdtd    total pressure in presence of the circular obstacle, 3d list [angle,distance,time].
    :param  p_f_fdtd    free-field pressure for the fdtd, 3d list [angle,distance,time].
    :param  f_fdtd      frequency sequence for the fdtd, 1d list of floats (Hz).
    :param  t_fdtd      time sequence for the fdtd, 1d list of floats (s).
    :param  Ts_fdtd     time step in the fdtd, float (s).
    :param  P_fdtd      Fourier transform of the scattered pressure for the FDTD, 3d list [angle,distance,frequency].
    :param  p_t_tlm     total pressure in presence of the circular obstacle, 3d list [angle,distance,time].
    :param  p_f_tlm     free-field pressure for the tlm, 3d list [angle,distance,time].
    :param  f_tlm       frequency sequence for the tlm, 1d list.
    :param  t_tlm       time sequence for the tlm, 1d list of floats (s).
    :param  Ts_tlm      time step in the tlm, float (s).
    :param  P_tlm       Fourier transform of the scattered pressure for the TLM, 3d list [angle,distance,frequency].

    :param  pan         analytic solution for plane waves scattered by a circular obstacle,
                        3d list [angle,distance,frequency].

    :param  Mag_Norm_fdtd   normalisation of the analytic solution magnitude in order to fit the numerical result
                            for the FDTD, float.
    :param  Mag_Norm_tlm    normalisation of the analytic solution magnitude in order to fit the numerical result
                            for the TLM, float.

    :param  two_norm_tlm    relative error in the 2-norm for the TLM as function
                            of frequency for each spatial step, 1d-array.
    :param  two_norm_fdtd   relative error in the 2-norm for the FDTD as function
                            of frequency for each spatial step, 1d-array.
    :param  max_norm_tlm    relative error in the max-norm for the TLM as function
                            of frequency for each spatial step, 1d-array.
    :param  max_norm_fdtd   relative error in the max-norm for the FDTD as function
                            of frequency for each spatial step, 1d-array.

    :param  ord_acc_tlm     order of accuracy between two consecutive grids in the 2-norm for the TLM
                            as function of frequency for each spatial step, 1d-array.
    :param  ord_acc_fdtd    order of accuracy between two consecutive grids in the 2-norm for the FDTD
                            as function of frequency for each spatial step, 1d-array.
    :param  ord_acc_tlm_max     order of accuracy between two consecutive grids in the max-norm for the TLM
                            as function of frequency for each spatial step, 1d-array.
    :param  ord_acc_fdtd_max    order of accuracy between two consecutive grids in the max-norm for the FDTD
                            as function of frequency for each spatial step, 1d-array.

    :param  freq_idx        chosen frequency index within the range(max_freq_idx), scalar.

    :return The log-log plot of the two-norm and max-norm, and the respective orders of accuracy
            for the chosen frequency.
    """
    f_idx_start = 1
    f_idx_end = 80  # 750Hz ; 120 # 1500 Hz

    for h_num, h in enumerate(h_set[:]):
        print 'h = %f ; %i/%i' % (h,h_num,len(h_set)-1)

        #   Load the numerical results and calculation of the FFTs
        res_path_fdtd = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case%i' % case, 'fdtd')
        p_f_fdtd = np.load(os.path.join(res_path_fdtd, 'p_h%i/p_%s.npy' % (h_num, 'f')))
        p_t_fdtd = np.load(os.path.join(res_path_fdtd, 'p_h%i/p_%s.npy' % (h_num, 't')))
        Ts_fdtd = np.load(os.path.join(res_path_fdtd, 'p_h%i/Ts.npy' % h_num))
        t_fdtd = np.load(os.path.join(res_path_fdtd, 'p_h%i/t.npy' % h_num))
        P_fdtd, f_fdtd = basic_fft(p_f_fdtd - p_t_fdtd, t_fdtd, Ts_fdtd, 2, 'hanning')
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'fdtd', 'p_h%i' % h_num)
        np.save(os.path.join(res_path, 'P_fdtd.npy'), P_fdtd)
        np.save(os.path.join(res_path, 'f_fdtd.npy'), f_fdtd)

        res_path_tlm = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case%i' % case, 'tlm')
        p_f_tlm = np.load(os.path.join(res_path_tlm, 'p_h%i/p_%s.npy' % (h_num, 'f')))
        p_t_tlm = np.load(os.path.join(res_path_tlm, 'p_h%i/p_%s.npy' % (h_num, 't')))
        Ts_tlm = np.load(os.path.join(res_path_tlm, 'p_h%i/Ts.npy' % h_num))
        t_tlm = np.load(os.path.join(res_path_tlm, 'p_h%i/t.npy' % h_num))
        P_tlm, f_tlm = basic_fft(p_f_tlm-p_t_tlm, t_tlm, Ts_tlm, 2, 'hanning')
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'tlm', 'p_h%i' % h_num)
        np.save(os.path.join(res_path, 'P_tlm.npy'), P_tlm)
        np.save(os.path.join(res_path, 'f_tlm.npy'), f_tlm)

        # Calculation of the analytic solutions
        dist_rcp_grid = np.load(os.path.join(res_path_fdtd, 'p_h%i/rcvdist.npy' % h_num))
        phi_rcp_grid = np.load(os.path.join(res_path_fdtd, 'p_h%i/rcvphi.npy' % h_num))
        ps1 = np.zeros((dist_rcp_grid.shape[0], len(dist_rcp_grid[0, :])), dtype=np.complex128)
        pan = np.zeros((dist_rcp_grid.shape[0], len(dist_rcp_grid[0, :]), len(f_fdtd[1:])), dtype=np.complex128)
        print 'frequency range: %.2f Hz - %.2f Hz.' %(f_fdtd[f_idx_start], f_fdtd[f_idx_end])

        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'analytic', 'p_h%i' % h_num)
        if not os.path.exists(res_path):
            os.makedirs(res_path)
            for i, f_value in enumerate(f_fdtd[f_idx_start:f_idx_end]):
                print 'Analytic calculation: takes time if a lot of frequencies!'
                Yb_ana = 1. * 10 ** -6 # rigid
                for m_ang in range(len(phi_rcp_grid[:, 0])):
                    for n_dist in range(len(dist_rcp_grid[0, :])):
                        ps1[m_ang, n_dist] = analytic_solution_scattered_pressure(rho, c, f_value, radius,
                                                                dist_rcp_grid[m_ang, n_dist],
                                                                (rho * c) / Yb_ana,
                                                                phi_rcp_grid[m_ang, n_dist])
                pan[:, :, i] = ps1
            np.save(os.path.join(res_path, 'p.npy'), pan)
