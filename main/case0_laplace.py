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

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fd_wave_eq')
site.addsitedir(fdtd_path)
from init_fd_modes import fd_helm_init

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc_fd_verif import error_calc


def main():
    """
    Two tests of the finite differences applied to the Helomholtz equation.
    1st test (case0): Laplace operator eigenfunction
    2nd test (case2): acoustic modes of a 2D rectangular plate
    :return: the plots of the convergence rate and order of accuracy
    """
    case = 0  # integer that sorts of the saved folders in the results dir.
    c = 340.  # sound speed, float (m.s-1).
    rho = 1.2  # air density, float (kg.m-3).
    fr = 50.
    wave_len = c / fr
    print wave_len / 10., wave_len / 12., wave_len / 180.

    # h_set = [wave_len / 10., wave_len / 20., wave_len / 60.]

    disp_inst_p = False  # display the instantaneous pressure, boolean.
    if case == 0 or case == 101:
        Lx = 1.
        Ly = 1.
        nx = 0
        ny = 0
        h_set = np.logspace(np.log10(0.001), np.log10(0.1), 50)

    if case == 2:
        Lx = 2. * 1.28
        Ly = 1.28
        nx = 1
        ny = 0
        h_set = np.logspace(np.log10(0.01), np.log10(0.16),
                            5)  # spatial step sequence, list of floats (m).
    for h_num, h in enumerate(h_set):
        fd_helm_init(h_num, h_set, Lx, Ly, nx, ny, fr, c, case,
                     disp_inst_p)
    error_calc(h_set, case)

if __name__ == '__main__':
    """
    dummy callback following the template
    """
    main()