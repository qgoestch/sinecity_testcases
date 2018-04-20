# -*- coding: utf-8 -*-
##
# \file     case1_geospr.py
# \title    Study of a Gaussian impulse geometrical spreading.
# \author   Pierre Chobeau
# \version  0.2
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 30 Aug.
##
import numpy as np
import os
import site


"""
.. module:: case1_geospr.py
   :platform: Unix, Windows
   :synopsis: Study of a Gaussian impulse geometrical spreading
   in a 2D domain using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- fdtd_srl_init_conv: initialization of the FDTD domain for the
study of geometrical spreading.

- tlm_srl_init_conv: initialization of the TLM domain for the
study of geometrical spreading.

- error_calc: results processing with FFT and errors calculations.

"""


base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_geospr import fdtd_srl_init_conv

tlm_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm_geospr import tlm_srl_init_conv

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc2_geospr import error_calc2
from errors_calc3_geospr import error_calc3
from spectral_error_geospr import spectral_error

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_errors_norms_cfa2018 import plot_errors_norms_cfa2018


def main(d_sr, f_max_src):
    """
    Each method (FDTD or TLM) is launched above a ground, then in free field.
    The numerical error is calculated in error_calc_geospr.py

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param f_max_src: approximated maximale frequency of the source signal (Hz).
    :type f_max_src: float
    """
    case = 1    # integer that sorts of the saved folders in the results dir.
    c = 340.    # sound speed, float (m.s-1).
    rho = 1.2   # air density, float (kg.m-3).
    T = 1./50 + 2 * 1./50. # T_delay + 2T_sim = simulation duration, float (s).
    # ~ T = 2.0 * 10 ** -2 s = 20 ms,
    # frequency res. df = 50 Hz, propagation distance = 6.8m
    # spatial step sequence, list of floats (m).
    h_set = np.logspace(np.log10(0.01), np.log10(0.16), 5)
    # time step sequence, list of floats (s).
    dt_set = np.logspace(np.log10(0.125), np.log10(2.0), 5) * 10**-4
    disp_inst_p = False     # display the instantaneous pressure, boolean.
    for h_idx, h in enumerate(h_set[:]):
        fdtd_srl_init_conv(h, h_idx, dt_set[h_idx], d_sr, T, f_max_src, rho, c,
                           case, disp_inst_p)
        tlm_srl_init_conv(h, h_idx, dt_set[h_idx], d_sr, T, f_max_src, rho, c,
                          case, disp_inst_p)

    # error_calc(d_sr, h_set[:], c, f_max_src, case)

    # error_calc2(d_sr, h_set[:], T, c, f_max_src, case)

    error_calc3(h_set, case)

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'tlm')
    one_norm_tlm = np.load(os.path.join(res_path, 'one_norm_tlm.npy'))
    two_norm_tlm = np.load(os.path.join(res_path, 'two_norm_tlm.npy'))
    max_norm_tlm = np.load(os.path.join(res_path, 'max_norm_tlm.npy'))
    ord_acc_one_tlm = np.load(os.path.join(res_path, 'ord_acc_one_tlm.npy'))
    ord_acc_two_tlm = np.load(os.path.join(res_path, 'ord_acc_two_tlm.npy'))
    ord_acc_max_tlm = np.load(os.path.join(res_path, 'ord_acc_max_tlm.npy'))

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'fdtd')
    one_norm_fdtd = np.load(os.path.join(res_path, 'one_norm_fdtd.npy'))
    two_norm_fdtd = np.load(os.path.join(res_path, 'two_norm_fdtd.npy'))
    max_norm_fdtd = np.load(os.path.join(res_path, 'max_norm_fdtd.npy'))
    ord_acc_one_fdtd = np.load(os.path.join(res_path, 'ord_acc_one_fdtd.npy'))
    ord_acc_two_fdtd = np.load(os.path.join(res_path, 'ord_acc_two_fdtd.npy'))
    ord_acc_max_fdtd = np.load(os.path.join(res_path, 'ord_acc_max_fdtd.npy'))

    plot_errors_norms_cfa2018(h_set, one_norm_tlm, one_norm_fdtd,
                      two_norm_tlm, two_norm_fdtd, max_norm_tlm, max_norm_fdtd,
                      ord_acc_one_tlm, ord_acc_one_fdtd,
                      ord_acc_two_tlm, ord_acc_two_fdtd,
                      ord_acc_max_tlm, ord_acc_max_fdtd,
                      case)

if __name__ == '__main__':
    """
    :param x_max: length of the receivers axis, float.
    :param x_min: first receiver location, float.
    :param d_sr: horizontal distance between the source and the receivers, 
    float (m).
    :param f_max_src: cutoff frequency of the Gaussian pulse, float (Hz).
    """
    x_max = 6.00
    x_min = 0.50
    d_sr = np.arange(x_min, x_max + .5, .5)
    print 'Number of points = %i' % len(d_sr)
    f_max_src = 2000.

    main(d_sr, f_max_src)
