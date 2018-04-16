# -*- coding: utf-8 -*-
##
# \file     case3_ground.py
# \title    Study of an acoustic impulse reflected by a ground.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 09 Aug.
##
import numpy as np
import os
import site


"""
.. module:: case3_ground.py
   :platform: Unix, Windows
   :synopsis: Study of an acoustic impulse reflected by a rigid ground in 
   a 2D domain using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- fdtd_srl_init_impgr: initialization of the FDTD domain for the
study of acoustic wave propagation above a reflecting ground.

- tlm_srl_init_impgr: initialization of the TLM domain for the
study of acoustic wave propagation above a reflecting ground.

- error_calc: results processing with FFT and errors calculations.

"""

fdtd_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_ground import fdtd_srl_init_impgr

tlm_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm_ground import tlm_srl_init_impgr

post_proc_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc_ground import error_calc


def main(d_sr, h_s, h_r, sigma, f_max_src):
    """
    Each method (FDTD or TLM) is launched above a ground, then in free field.
    The numerical error is calculated in error_calc_ground.py

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param sigma: pecific airflow resistivity (kNm-4s==CGS).
    :type sigma: float
    :param f_max_src: approximated maximale frequency of the source signal -
    Gaussian pulse - (Hz).
    :type f_max_src: float

    :param case: integer that sorts of the saved folders in the results dir.
    :param c: sound speed, float (m.s-1).
    :param rho: air density, float (kg.m-3).
    :param T: simulation duration, float (s).
    :param h_set: spatial step sequence, list of floats (m).
    :param dt_set: time step sequence, list of floats (s).
    :param disp_inst_p: display the instantaneous pressure, boolean.
    """
    case = 3
    c = 340.
    rho = 1.2
    T = 1. / 25.  # 0.04 s --> 13.6 m propag. dist.
    # (set T_delay in init_*.py accordingly: T/2) ; df=25Hz
    # T = 1. / 50.  # 0.02 s --> 6.8 m propag. dist.
    # (set T_delay in init_*.py accordingly: T/4) ; df=50Hz

    h_set = np.logspace(np.log10(0.01), np.log10(0.16), 5)
    dt_set = np.logspace(np.log10(0.125), np.log10(2.0), 5)*10**-4
    disp_inst_p = False

    # for h_idx, h in enumerate(h_set[:]):
    #     for ff in [False, True]:
    #         fdtd_srl_init_impgr(dt_set[h_idx], h, h_idx, h_set, d_sr, h_s, h_r,
    #                             T, f_max_src, rho, sigma, case, ff, disp_inst_p)
    #         tlm_srl_init_impgr(dt_set[h_idx], h, h_idx, h_set, d_sr, h_s, h_r,
    #                            T, f_max_src, rho, sigma, case, ff, disp_inst_p)

    error_calc(d_sr, h_s, h_r, h_set[:], rho, c, sigma, case,
               disp_att_spect=False, disp_errors=True)

if __name__ == '__main__':
    two_coarsest_spatial_step = 2 * 0.16
    h_s = 2. * two_coarsest_spatial_step
    x_max = 16. * two_coarsest_spatial_step
    x_min = two_coarsest_spatial_step
    y_max = 7. * two_coarsest_spatial_step
    y_min = two_coarsest_spatial_step
    d_sr = np.arange(x_min, x_max + two_coarsest_spatial_step,
                     two_coarsest_spatial_step)
    h_r = np.arange(y_min, y_max + two_coarsest_spatial_step,
                    two_coarsest_spatial_step)
    sigma = 20000
    f_max_src = 2000.

    main(d_sr, h_s, h_r, sigma, f_max_src)
