# -*- coding: utf-8 -*-
##
# \file     init_fd_modes.py
# \title    Set the 2D domain for FD method applied on the wave equation.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2018, 30 Jan.
##
import numpy as np
import os
import site
import scipy.special as sp

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fd_core')
site.addsitedir(fdtd_core_path)
from upd_fd import upd_p_fd_helm_bc_val_pb, upd_p_fd_helm_sommerfeld_bc, \
                    upd_p_fd_helm_sommerfeld_bc_source

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_modes, \
    green_sol_to_2D_helmholtz_eq, testcase1_eigenfunct, \
    analytic_solution_modes_fd_helmholtz

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure_fd_method
from plot_errors_norms import plot_error_basic

post_proc_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'post_proc')
site.addsitedir(post_proc_path)
from errors_calc_fd_verif import error_calc


def fd_helm_init(h_num, h_set, Lx, Ly, nx, ny, fr, c, case, disp_inst_p):
    """
    Setting the 2D geometries and running the FD method on the Helmholtz
    equation.

    :param h_num: spatial step index.
    :type h_num: int
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param Lx: length of the box followinf the x axis (m)
    :type Lx: float
    :param Ly: length of the box followinf the y axis (m)
    :type Ly: float
    :param nx: mode number following the x direction
    :type nx: int
    :param ny: mode number following the y direction
    :type ny: int
    :param fr: frequency of the signal (Hz).
    :type fr: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param disp_inst_p: display the instantaneous pressure.
    :type disp_inst_p: bool
    :return: the acoustic pressure at the pre-defined receivers' locations.
    :rtype: (1+1)D array of floats 64
    """
    # =========================================================================
    #   Grid parameters
    # =========================================================================
    k = 2. * np.pi * fr / c       # wave number (rad.m-1)
    Lx = Lx + h_set[h_num]
    Ly = Ly + h_set[h_num]
    Nx = np.int(np.round(Lx / h_set[h_num]))  # number of node following x axis
    Ny = np.int(np.round(Ly / h_set[h_num]))  # number of node following y axis
    x = np.linspace(0, Lx, Nx + 1)  # discretized x axis
    y = np.linspace(0, Ly, Ny + 1)  # discretized y axis
    dx = np.float64(x[1] - x[0])    # spatial step for the x direction
    dx = round(dx, 5)
    dy = dx

    print '                 FD for Helmholtz in 2D box                    '
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: h=%g m.' % dx
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' \
          % (Nx, Ny, Lx, Ly)

    # =========================================================================
    #   Variables
    # =========================================================================
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p_an = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    f_source = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)

    # ========================================================================
    #   Definition of the initial conditions: centered point source
    # ========================================================================
    # xs_idx = int(p.shape[0] / 2.)
    # ys_idx = int(p.shape[1] / 2.)
    # f_source[:, :] = - green_sol_to_2D_helmholtz_eq(Nx, Ny, xs_
#                                                     idx, ys_idx, dx, k)

    # =========================================================================
    #   Test case 0 from [Eq. (39), sutmann_jcam2007]
    # =========================================================================
    if case == 0:
        xv = x[:, np.newaxis] - (dx / 2)
        yv = y[np.newaxis, :] - (dx / 2)
        f_source[1:-1, 1:-1] = (3. * np.pi**2 / Lx + 1j * k**2) * \
                                np.sin(np.pi * xv[1:-1, :] / Lx) * \
                                np.sin(np.pi * yv[:, 1:-1] / Ly)
        p_an[1:-1, 1:-1] = testcase1_eigenfunct(Lx, Ly, dx, x, y)

    # =========================================================================
    #   Test case 3 from [Eq. (40), sutmann_jcam2007] or [manohar_jcomphy1983]
    # =========================================================================
    if case == 101:
        xv = x[:, np.newaxis] - (dx / 2)
        yv = y[np.newaxis, :] - (dx / 2)
        p[1:-1, 1:-1] = 1. / np.cosh(10.) * (
                        np.cosh(10. * xv[1:-1, :]) +
                        np.cosh(10. * yv[:, 1:-1]))

    # =========================================================================
    #   Case 4: acoustic MODES in a 2D box
    # =========================================================================
    if case == 2:
        f_source[1:-1, 1:-1] = analytic_solution_modes_fd_helmholtz(dx, nx, ny,
                                                                    Lx, Ly,
                                                                    Nx, Ny,
                                                                    x, y)
        p_an[1:-1, 1:-1] = analytic_solution_modes_fd_helmholtz(dx, nx, ny,
                                                                Lx, Ly,
                                                                Nx, Ny, x, y)

    # =========================================================================
    #   Calculation of the pressure : solver
    # =========================================================================
    # tol = 9.9999*10**-1  # tolerance parameters for establishing the pressure field
    # i = 0  # loop count
    # max_v_diff = 1
    # while max_v_diff > tol:
    #     i += 1
    #     print i
    #     p_old = p.copy()
    #     p = upd_p_fd_helm_sommerfeld_bc_source(p, f_source, k, dx)
    #     # p = upd_p_fd_helm_bc_val_pb(p, k, dx)
    #
    #     # check for convergence
    #     p_diff = p - p_old
    #     max_v_diff = np.absolute(p_diff).max()

    p = upd_p_fd_helm_sommerfeld_bc_source(p, f_source, k, dx)

    p_normalized_fd = np.real(p) / max(abs(np.max(np.real(p[1:-1, 1:-1]))),
                                       abs(np.min(np.real(p[1:-1, 1:-1]))))
    p_normalized_an = np.real(p_an) / max(abs(np.max(np.real(p_an[1:-1, 1:-1]))),
                                          abs(np.min(np.real(p_an[1:-1, 1:-1]))))

    if disp_inst_p:
        # instatenous_pressure_fd_method(p_normalized_an, dx, Lx, Ly, case, False)
        instatenous_pressure_fd_method(p_normalized_fd, dx, Lx, Ly, case, False)

    print np.max(np.real(p[1:-1, 1:-1])), np.min(np.real(p[1:-1, 1:-1]))

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i'
                            % case, 'fd')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    np.save(os.path.join(res_path, 'p_fd_%i.npy' % h_num), p_normalized_fd)
    np.save(os.path.join(res_path, 'p_an_%i.npy' % h_num), p_normalized_an)


if __name__ == '__main__':
    case = 0    # integer that sorts of the saved folders in the results dir.
    c = 340.    # sound speed, float (m.s-1).
    rho = 1.2   # air density, float (kg.m-3).
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
        Lx = 2.*1.28
        Ly = 1.28
        nx = 1
        ny = 0
        h_set = np.logspace(np.log10(0.01), np.log10(0.16),
                            5)  # spatial step sequence, list of floats (m).
    for h_num, h in enumerate(h_set):
        fd_helm_init(h_num, h_set, Lx, Ly, nx, ny, fr, c, case,
                                disp_inst_p)
    error_calc(h_set, case)
