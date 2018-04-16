# -*- coding: utf-8 -*-
##
# \file     paumond_facade.py
# \title    Study of an acoustic impulse in half street (sidewalk+facade).
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 15 Jan.
##
import numpy as np
import os
import site


"""
.. module:: paumond_facade.py
   :platform: Unix, Windows
   :synopsis: Study of an acoustic impulse in half street (sidewalk+facade) in 
   a 2D domain using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- init_fdtd_paumond_facade: initialization of the FDTD domain for the
study of acoustic wave propagation above a reflecting ground.

- init_tlm_paumond_facade: initialization of the TLM domain for the
study of acoustic wave propagation above a reflecting ground.

- level_calc_paumond_facade_2: results processing with FFT and 1/3rd octave
                                box plots.

"""

fdtd_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_paumond_facade import fdtd_srl_init

fdtd_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(fdtd_path)
from init_tlm_paumond_facade import tlm_srl_init

post_proc_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from level_calc_paumond_facade_2 import levels


def main(d_sr, d_sr_ref, h_s, h_r, f_max_src):
    """
    Each method (FDTD or TLM) is launched. The numerical error is calculated
    in level_calc_paumond_facade_2.py, and presented in boxplots.

    :param d_sr: horizontal distances between the facade and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
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
    case = 668
    c = 340.
    rho = 1.2
    T_total = 1. / 10.  # 0.01 s --> propa. dist.: 34.0 m
    T_delay = 0.02
    disp_inst_p = False

    fdtd_srl_init(d_sr, d_sr_ref, h_s, h_r, T_total, T_delay, f_max_src,
                  rho, c, case, disp_inst_p)
    tlm_srl_init(d_sr, d_sr_ref, h_s, h_r, T_total, T_delay, f_max_src,
                  rho, c, case, disp_inst_p)

    levels(d_sr, d_sr_ref, h_r, case, 'fdtd')

from random import gauss
if __name__ == '__main__':
    d_sr = np.arange(0.10, 0.30, 0.01)  # dist. facade - mic. balcons

    # d_sr_ref = np.linspace(1.5, 3.0, len(d_sr))   # dist. facade - mic. piéton
    d_sr_ref = [gauss(2.25, 0.75) for i in range(len(d_sr))]

    print len(d_sr), len(d_sr_ref)

    # h_r = [1.5, 4.0, 7.0, 12.0]
    h_r = [1.5, 3., 4., 4., 4., 5., 7., 7., 7., 7., 7., 9., 9., 9., 9.,
           12., 12., 12., 15., 15., 15., 15.]  # hauteur des balcons
    h_s = 0.20
    f_max_src = 2244. # 1.122 * 2000. to get the 1/3rd octave 2000 Hz

    main(d_sr, d_sr_ref, h_s, h_r, f_max_src)
