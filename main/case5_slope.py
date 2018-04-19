# -*- coding: utf-8 -*-
##
# \file     case5_slope.py
# \title    Study of an impulse diffracted by a corner from an upward sloping part.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 27 Sep.
##
import numpy as np
import os
import site


"""
.. module:: case5_slope.py
   :platform: Unix, Windows
   :synopsis: Study of an acoustic impulse diffracted by a corner from 
   an upward sloping part in a 2D domain using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- fdtd_srl_init_slope: initialization of the FDTD domain for the
study of acoustic wave propagation above a corner.

- tlm_srl_init_slope: initialization of the TLM domain for the
study of acoustic wave propagation above a corner.

- error_calc: results processing with FFT and errors calculations.

"""


base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_upslope import fdtd_srl_init_slope

tlm_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm_upslope import tlm_srl_init_slope

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc_slope import error_calc


def main(freq, h_s, l1max, l2max, zmax, alpha):
    """
    Each method (FDTD or TLM) is launched above a ground, then in free field.
    The numerical error is calculated in error_calc_ground.py
    :param freq: frequency of the source signal (Hz).
    :type freq: float
    :param h_s: source heights (m).
    :type h_s: float
    :param l1max: length of the horizontal part (m).
    :type l1max: float
    :param l2max: length of the upward-sloping part (m).
    :type l2max: float
    :param zmax: height of the the two parts (m).
    :type zmax: float
    :param alpha: angle of the upward-sloping part (rad).
    :type alpha: float

    :param case: integer that sorts of the saved folders in the results directory.
    :param c: sound speed, float (m.s-1).
    :param rho: air density, float (kg.m-3).
    :param T: simulation duration, float (s).
    :param h_set: spatial step sequence, list of floats (m).
    :param disp_inst_p: display the instantaneous pressure, boolean.
    """
    case = 5
    c = 340.
    rho = 1.2
    # T = 1. / 25.  # 0.04 seconds, df=25Hz, propag. dist. = 13.6m
    T = 4. / 10.  # 0.4 seconds, df=2.5Hz, propag. dist. = 136m
    # h_set = np.logspace(np.log10(0.01), np.log10(0.16), 5)
    h_set = [0.10, 0.09]
    # disp_inst_p = True
    # for al in alpha:
    #     for h_idx, h in enumerate(h_set[:]):
    #         fdtd_srl_init_slope(freq, h_s, l1max, l2max, zmax, al, c,
    #                             h, h_idx + 4, h_set, T, rho, case, disp_inst_p)
    #         tlm_srl_init_slope(freq, h_s, l1max, l2max, zmax, al, c,
    #                             h, h_idx + 4, h_set, T, rho, case, disp_inst_p)

    error_calc(h_set, h_s, alpha, case, 'tlm')

if __name__ == '__main__':
    freq = 800.
    h_s = 2.
    l1max = 60.
    l2max = 40.
    zmax = 40.
    alpha = [10., 20.]
    main(freq, h_s, l1max, l2max, zmax, alpha)