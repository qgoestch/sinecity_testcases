# -*- coding: utf-8 -*-
##
# \file     case2_modes.py
# \title    Study of a 2D domain acoustic modes.
# \author   Pierre Chobeau
# \version  0.2
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 25 Jul.
##
import numpy as np
import os
import site


"""
.. module:: case2_modes.py
   :platform: Unix, Windows
   :synopsis: Study of the acoustic modes in a 2D domain 
   made of reflecting boundary conditions using the FDTD and the TLM methods.

.. moduleauthor:: Pierre Chobeau <pierre.chobeau@ifsttar.fr>

List of required functions
==========================

- fdtd_srl_init_conv_box: initialization of the FDTD domain for the
study of the acoustic modes using initial conditions.

- tlm_srl_init_conv_box: initialization of the TLM domain for the
study of the acoustic modes using initial conditions.

- error_calc: results processing with FFT and errors calculations.

"""


base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd_modes import fdtd_srl_init_conv_box

tlm_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm_modes import tlm_srl_init_conv_box

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc2_modes import error_calc2

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_errors_norms_cfa2018 import plot_errors_norms_cfa2018, plot_errors_norms_fd_fdtd_tlm_cfa2018


def main(Lx, Ly, nx, ny, fr):
    """
    Each method (FDTD or TLM) is called using the init_*.py functions.
    The resutls are post-processed using error_calc function.
    
    :param Lx: length of the box followinf the x axis (m)
    :type Lx: float
    :param Ly: length of the box followinf the y axis (m)
    :type Ly: float
    :param nx: mode number following the x direction
    :type nx: int
    :param ny: mode number following the y direction
    :type ny: int
    :param fr: frequency of the sine source signal (Hz)
    :type fr: float

    """
    case = 2    # integer that sorts of the saved folders in the results directory.
    c = 340.    # sound speed, float (m.s-1).
    rho = 1.2   # air density, float (kg.m-3).
    dt_coarse = 2. * 10 ** -4   # corsest time step that corresponds to the coarsest grid, float (s).
    n_it_save = 2   # number of iteration before saving the numerical pressure relative to the coarsest grid, int.
    t_save = n_it_save * dt_coarse  # exact time at which the pressures are compared, int.
    T = (n_it_save + 2) * dt_coarse # simulation duration, float (s).
    h_set = np.logspace(np.log10(0.01), np.log10(0.16), 5)  # spatial step sequence, list of floats (m).
    disp_inst_p = False  # display the instantaneous pressure, boolean.

    for h_num,h in enumerate(h_set):
        fdtd_srl_init_conv_box(h, h_num, h_set, t_save, T, Lx, Ly, nx, ny, fr,
                               rho, c, case, disp_inst_p)
        tlm_srl_init_conv_box(h, h_num, h_set, t_save, T, Lx, Ly, nx, ny, fr,
                              rho, c, case, disp_inst_p)

    # error_calc(h_set, case)

    error_calc2(h_set, case)

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

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'fd')
    one_norm_fd = np.load(os.path.join(res_path, 'one_norm_fd.npy'))
    two_norm_fd = np.load(os.path.join(res_path, 'two_norm_fd.npy'))
    max_norm_fd = np.load(os.path.join(res_path, 'max_norm_fd.npy'))
    ord_acc_one_fd = np.load(os.path.join(res_path, 'ord_acc_one_fd.npy'))
    ord_acc_two_fd = np.load(os.path.join(res_path, 'ord_acc_two_fd.npy'))
    ord_acc_max_fd = np.load(os.path.join(res_path, 'ord_acc_max_fd.npy'))

    # plot_errors_norms(h_set, one_norm_tlm, one_norm_fdtd,
    #                   two_norm_tlm, two_norm_fdtd, max_norm_tlm, max_norm_fdtd,
    #                   ord_acc_one_tlm, ord_acc_one_fdtd,
    #                   ord_acc_two_tlm, ord_acc_two_fdtd,
    #                   ord_acc_max_tlm, ord_acc_max_fdtd,
    #                   case)

    plot_errors_norms_fd_fdtd_tlm_cfa2018(h_set, one_norm_fd, one_norm_tlm, one_norm_fdtd,
                                  two_norm_fd, two_norm_tlm, two_norm_fdtd,
                                  max_norm_fd, max_norm_tlm, max_norm_fdtd,
                                  ord_acc_one_fd, ord_acc_one_tlm, ord_acc_one_fdtd,
                                  ord_acc_two_fd, ord_acc_two_tlm, ord_acc_two_fdtd,
                                  ord_acc_max_fd, ord_acc_max_tlm, ord_acc_max_fdtd, case)

if __name__ == '__main__':
    Lx = 2.*1.28
    Ly = 1.28
    fr = 358.391
    nx = 2
    ny = 1
    main(Lx, Ly, nx, ny, fr)