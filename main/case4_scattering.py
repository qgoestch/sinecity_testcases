# -*- coding: utf-8 -*-
##
# \file     case4_scattering.py
# \title    Study of an acoustic impulse scattered by a circular obstacle.
# \author   Pierre Chobeau
# \version  0.2
# \date     2017, 07 Sep.
##
import numpy as np
import os
import site


"""
.. module:: case4_scattering.py
   :platform: Unix, Windows
   :synopsis: Study of an acoustic impulse scattered by a circular obstacle in 
   a 2D domain using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- fdtd_srl_init_scat: initialization of the FDTD domain for the
study of acoustic wave scattered by a circular obstacle.

- tlm_srl_init_scat: initialization of the TLM domain for the
study of acoustic wave scattered by a circular obstacle.

- error_calc: results processing with FFT and errors calculations.

"""


base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_scat import fdtd_srl_init_scat

tlm_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm_scat import tlm_srl_init_scat

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc2_scat import error_calc2


def main(f_max_src):
    """
    :param f_max_src: approximated maximale frequency of the source signal -
    Gaussian pulse - (Hz).
    :type f_max_src: float
    
    :param case: integer that sorts of the saved folders in the results dir.
    :param c: sound speed, float (m.s-1).
    :param rho: air density, float (kg.m-3).
    :param T_sim: simulation duration after the pulse starts, float (s).
    :param T_delay: time before the pulse starts - zero pressure imposed, float (s).
    :param T: simulation duration, float (s).
    :param radius: radius of the scatterer (m)
    :param h_set: spatial step sequence, list of floats (m).
    :param disp_inst_p: display the instantaneous pressure, boolean.
    """
    case = 4
    c = 340.
    rho = 1.2
    # T = 1. / 50.    # 0.02 seconds, df=50Hz, propag. dist. = 6.8m
    # T = 1. / 25.  # 0.04 seconds, df=25Hz, propag. dist. = 13.6m
    T_sim  = 0.045
    T_delay= 0.035
    T      = T_sim + T_delay
    
    radius = 0.301
    h_set = [0.0213, 0.0251, 0.0274, 0.0355, 0.0405, 0.0430,
             0.0475, 0.0495, 0.0550, 0.0580, 0.0670, 0.0695]
    disp_inst_p = False

    # for h_idx, h in enumerate(h_set[:]):
    #     for ff in [False, True]:
    #         fdtd_srl_init_scat(h, h_idx, radius, T_delay, T, f_max_src, rho, c,
    #                            case, disp_inst_p,ff)
    #         tlm_srl_init_scat(h, h_idx, radius, T_delay, T, f_max_src, rho, c,
    #                           case, disp_inst_p,ff)

    # error_calc(h_set, rho, c, radius, case)
    error_calc2(h_set, rho, c, radius, case)


if __name__ == '__main__':
    f_max_src = 2000.
    main(f_max_src)
