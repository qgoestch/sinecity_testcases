# -*- coding: utf-8 -*-
##
# \file     upd_fd.py
# \title    FD updates of the Helmholtz equation - frequency domain.
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 01 Feb.
##

import numpy as np


def upd_p_fd_helm_sommerfeld_bc_source(p, f, k, dx):
    """
    Pressure update of the 2D Helmholtz equation with Sommerfeld boundary
    condition, see [Eqs. (4), (5), hegedus_atj2010].

    :param p: updated pressure at time n+1 (Pa).
    :type p: 2D numpy array of floats.
    :param f: source term (Pa).
    :type f: 2D numpy array of floats.
    :param k: wavenumber (rad.m-1).
    :type k: float
    :param dx: spatial step (m).
    :type dx: float

    :return: the pressure from the discretized helmholtz eq., variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   main solver (inside the domain)
    # =========================================================================
    p[1:-1, 1:-1] = 1. / ((k ** 2 * dx ** 2) - 4.) * \
                    (- f[1:-1, 1:-1] -
                    (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]))

    # =========================================================================
    #   B.C. on the frame Sommerfeld
    # =========================================================================
    # edges
    for i in range(1, p.shape[0]-2):
        j_min = 1
        p[i, j_min] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i, j_min + 1] + p[i-1, j_min] + p[i+1, j_min])
        j_max = p.shape[1]-2
        p[i, j_max] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i, j_max - 1] + p[i-1, j_max] + p[i+1, j_max])

    for j in range(1, p.shape[1]-2):
        i_min = 1
        p[i_min, j] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i_min + 1, j] + p[i_min, j-1] + p[i_min, j+1])
        i_max = p.shape[0]-2
        p[i_max, j] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i_max - 1, j] + p[i_max, j-1] + p[i_max, j+1])

    # corners
    p[1, 1] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[1, 2] + 2 * p[2, 1])
    p[1, p.shape[1]-2] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[1, p.shape[1]-3] + 2 * p[2, p.shape[1]-2])
    p[p.shape[0]-2, 1] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[p.shape[0]-3, 1] + 2 * p[p.shape[0]-2, 2])
    p[p.shape[0]-2, p.shape[1]-2] = \
        (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
        (2. * p[p.shape[0]-3, p.shape[1]-2] + 2 * p[p.shape[0]-2, p.shape[1]-3])
    return p


def upd_p_fd_helm_sommerfeld_bc(p, k, dx):
    """
    Pressure update of the 2D Helmholtz equation with Sommerfeld boundary
    condition, see [Eqs. (4), (5), hegedus_atj2010].

    :param p: updated pressure at time n+1, (Pa).
    :type p: 2D numpy array of floats.
    :param k: wavenumber (rad.m-1).
    :type k: float
    :param dx: spatial step (m).
    :type dx: float

    :return: the pressure from the discretized helmholtz eq., variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   main solver (inside the domain)
    # =========================================================================
    p[1:-1, 1:-1] = 1. / (4. - (k ** 2 * dx ** 2)) * \
                    (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2])

    # =========================================================================
    #   B.C. on the frame Sommerfeld
    # =========================================================================
    # edges
    for i in range(1, p.shape[0]-2):
        j_min = 1
        p[i, j_min] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i, j_min + 1] + p[i-1, j_min] + p[i+1, j_min])
        j_max = p.shape[1]-2
        p[i, j_max] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i, j_max - 1] + p[i-1, j_max] + p[i+1, j_max])

    for j in range(1, p.shape[1]-2):
        i_min = 1
        p[i_min, j] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i_min + 1, j] + p[i_min, j-1] + p[i_min, j+1])
        i_max = p.shape[0]-2
        p[i_max, j] = (1. / (4. - 2.*1j*k*dx - k**2 * dx**2)) * \
                      (2. * p[i_max - 1, j] + p[i_max, j-1] + p[i_max, j+1])

    # corners
    p[1, 1] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[1, 2] + 2 * p[2, 1])
    p[1, p.shape[1]-2] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[1, p.shape[1]-3] + 2 * p[2, p.shape[1]-2])
    p[p.shape[0]-2, 1] = (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
                  (2. * p[p.shape[0]-3, 1] + 2 * p[p.shape[0]-2, 2])
    p[p.shape[0]-2, p.shape[1]-2] = \
        (1. / (4. - 2. * 1j * k * dx - k ** 2 * dx ** 2)) * \
        (2. * p[p.shape[0]-3, p.shape[1]-2] + 2 * p[p.shape[0]-2, p.shape[1]-3])
    return p


def upd_p_fd_helm_bc_val_pb(p, k, dx):
    """
    Finite difference for the Helmholtz equation solving two different
    boundary value problems found in [zhang_fd_method_for_helmholtz].

    :param p: updated pressure at time n+1, (Pa).
    :type p: 2D numpy array of floats.
    :param k: wavenumber (rad.m-1).
    :type k: float
    :param dx: spatial step (m).
    :type dx: float

    :return: the pressure from the discretized helmholtz eq., variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   main solver (inside the domain)
    # =========================================================================
    p[1:-1, 1:-1] = 1. / (4. - (k ** 2 * dx ** 2)) * \
                    (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2])

    # =======================================================================
    #   B.C. see [Eq. (2.22), zhang_fd_method]
    #   Seems ok for : 100 Hz, lambda/40, tol = 5.35*1e-3
    # =======================================================================
    for i_idx, i_val in enumerate(range(-int((p.shape[0] - 1) / 2),
                                        int(p.shape[0] / 2))):
        p[i_idx, 0] = np.sin(np.pi * i_val * dx / 6.)
        p[i_idx, p.shape[1]-1] = np.sin(np.pi * i_val * dx / 6.)
    for j_idx, j_val in enumerate(range(-int((p.shape[1] - 1) / 2),
                                        int(p.shape[1] / 2))):
        p[0, j_idx] = np.sin(np.pi * -int((p.shape[0] - 1) / 2) * dx / 6.)
        p[p.shape[0]-1, j_idx] = np.sin(np.pi * int(p.shape[0] / 2) * dx / 6.)

    # # =======================================================================
    # #   B.C. see [Eq. (2.23), zhang_fd_method]
    # #   Seems ok for : 100 Hz, lambda/40, tol = 0.47*1e-3
    # # =======================================================================
    # for i_idx, i_val in enumerate(range(-int((p.shape[0] - 1) / 2),
    #                                     int(p.shape[0] / 2) + 1)):
    #     p[i_idx, 0] = -1. / ((i_val * dx) ** 2 +
    #                          (-int(p.shape[1] / 2) * dx) ** 2)
    #     p[i_idx, p.shape[1]-1] = -1. / ((i_val * dx) ** 2 +
    #                                     (int(p.shape[1] / 2) * dx) ** 2)
    # for j_idx, j_val in enumerate(range(-int((p.shape[1] - 1) / 2),
    #                                     int(p.shape[1] / 2) + 1)):
    #     p[0, j_idx] = -1. / ((-int(p.shape[0] / 2) * dx) ** 2 +
    #                          (j_val * dx) ** 2)
    #     p[p.shape[0]-1, j_idx] = -1. / ((int(p.shape[0] / 2) * dx) ** 2 +
    #                                     (j_val * dx) ** 2)

    # # =======================================================================
    # #   B.C. on the frame Dirichlet = perfectly reflecting.
    # # =======================================================================
    # j = 1;              p[:, j] = p[:, j+1]
    # j = p.shape[1]-2;   p[:, j] = p[:, j-1]
    # i = 1;              p[i, :] = p[i+1, :]
    # i = p.shape[0]-2;   p[i, :] = p[i-1, :]


    return p


def upd_p_fdtd_srl_2D_scat(p, p1, p2, fsrc, fsrc2, Nb, c, rho, Ts, dx, Cn,
                           A, B, C, x_in_idx, y_in_idx,
                           x_edges_idx, y_edges_idx,
                           x_corners_idx, y_corners_idx):
    """
    This FDTD update is designed for case 4: scattering.
    It calculates the pressure at the discrete time n+1 (p),
    as a function of the pressures at time n (p1) and n-1 (p2).
    The stencil follows a Standard RectiLinear (SRL) implementation.
    This update is implemented with a frequency-dependent boundary condition.

    :param p: updated pressure at time n+1, (Pa).
    :type p: 2D numpy array of floats.
    :param p1: current pressure at time n, (Pa).
    :type p1: 2D numpy array of floats.
    :param p2: past pressure at time n-1, (Pa).
    :type p2: 2D numpy array of floats.
    :param fsrc: soft source at time n+1, (Pa).
    :type fsrc: 2D numpy array of floats.
    :param fsrc1: soft source at time n, (Pa).
    :type fsrc1: 2D numpy array of floats.
    :param fsrc2: soft source at time n-1, (Pa).
    :type fsrc2: 2D numpy array of floats.
    :param Nb: boundary of the domain (1 if BC, 0 else).
    :type Nb: 2D numpy array of integers.
    :param c: sound speed (m.s-1).
    :type c: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param Ts: time step (s).
    :type Ts: float
    :param dx: spatial step (m).
    :type dx: float
    :param Cn: Courant number.
    :type Cn: float
    :param A: inertance of the boundary.
    :type A: float
    :param B: stiffness of the boundary.
    :type B: float
    :param C: resistivity of the boundary.
    :type C: float
    :param x_in_idx: x coordinates of the p=0 cells under the slope BC cells
    :type x_in_idx: list of integers
    :param y_in_idx: y coordinates of the p=0 cells under the slope BC cells
    :type y_in_idx: list of integers
    :param x_edges_idx: x coordinates of the slope BC cells that have single
    face contact, i.e. edges
    :type x_edges_idx: list of integers
    :param y_edges_idx: y coordinates of the slope BC cells that have single
    face contact, i.e. edges
    :type y_edges_idx: list of integers
    :param x_corners_idx: x coordinates of the slope BC cells that have two
    faces contact, i.e. corners
    :type x_corners_idx: list of integers
    :param y_corners_idx: y coordinates of the slope BC cells that have two
    faces contact, i.e. corners
    :type y_corners_idx: list of integers
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    j = 0;
    p[:, j] = 0
    j = p.shape[1] - 1;
    p[:, j] = 0
    i = 0;
    p[i, :] = 0
    i = p.shape[0] - 1;
    p[i, :] = 0
    # =========================================================================
    #   main solver (inside the domain): SRL FDTD <==> Cubic cell FVTD
    # =========================================================================
    A_corr_src = -1. / 3.1739566e3
    p[1:-1, 1:-1] = \
    1. / (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) + (.5 * B)))) * \
    (
        (
            (Cn ** 2) *
            (p1[2:, 1:-1] + p1[:-2, 1:-1] +
             p1[1:-1, 2:] + p1[1:-1, :-2] -
             ((4. - Nb[1:-1, 1:-1]) * p1[1:-1, 1:-1]))
        )
        + (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
        p1[1:-1, 1:-1]
        - (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) - (.5 * B)))) *
        p2[1:-1, 1:-1]
    ) - A_corr_src * rho * (c ** 2 * Ts ** 2 / dx ** 2) * \
        (fsrc[1:-1, 1:-1] - fsrc2[1:-1, 1:-1]) / (2 * Ts)

    # =========================================================================
    #   Zero pressure for the cells inside the cylinder
    # =========================================================================
    for ndx_0 in range(len(x_in_idx)):
        p[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0

    # =========================================================================
    #   Hamilton's update in DAFx 2014 Eq. (39) with "full cell" implementation
    # =========================================================================
    #   Edge-cell
    for ndx_ed in range(len(x_edges_idx)):
        i = x_edges_idx[ndx_ed]
        j = y_edges_idx[ndx_ed]
        p[i, j] = \
        1. / (1. + (Nb[i, j] * Cn * ((A / Ts) + (.5 * B)))) * \
        (   (   (Cn ** 2) *
                (p1[i + 1, j] + p1[i - 1, j] +
                 p1[i, j + 1] + p1[i, j - 1] -
                 (3.* p1[i, j])   )       )
            + (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
            p1[i, j]
            - (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * B)))) *
            p2[i, j]
        )

    #   Corner-cell
    for ndx_co in range(len(x_corners_idx)):
        i = x_corners_idx[ndx_co]
        j = y_corners_idx[ndx_co]
        p[i, j] = \
        1. / (1. + (Nb[i, j] * Cn * ((A / Ts) + (.5 * B)))) * \
        (   (   (Cn ** 2) *
                (p1[i + 1, j] + p1[i - 1, j] +
                 p1[i, j + 1] + p1[i, j - 1] -
                 (2.* p1[i, j])   )       )
            + (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
            p1[i, j]
            - (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * B)))) *
            p2[i, j]
        )

    # =========================================================================
    #   B.C. on the frame Dirichlet = perfectly reflecting.
    # =========================================================================
    j = 1;
    p[:, j] = p[:, j + 1]
    j = p.shape[1] - 2;
    p[:, j] = p[:, j - 1]
    i = 1;
    p[i, :] = p[i + 1, :]
    i = p.shape[0] - 2;
    p[i, :] = p[i - 1, :]
    return p


def upd_p_fdtd_srl_2D_slope(p, p1, p2, fsrc, fsrc2, Nb, c, rho, Ts, dx, Cn,
                            A, B, C, x_in_idx, y_in_idx, x_edges_idx,
                            y_edges_idx, x_corners_idx,
                            y_corners_idx, slope_start):
    """
     This FDTD update is designed for case 5: slope.
    It calculates the pressure at the discrete time n+1 (p),
    as a function of the pressures at time n (p1) and n-1 (p2).
    The stencil follows a Standard RectiLinear (SRL) implementation.
    This update is implemented with a frequency-dependent boundary condition.

    :param p: updated pressure at time n+1, (Pa).
    :type p: 2D numpy array of floats.
    :param p1: current pressure at time n, (Pa).
    :type p1: 2D numpy array of floats.
    :param p2: past pressure at time n-1, (Pa).
    :type p2: 2D numpy array of floats.
    :param fsrc: soft source at time n+1, (Pa).
    :type fsrc: 2D numpy array of floats.
    :param fsrc1: soft source at time n, (Pa).
    :type fsrc1: 2D numpy array of floats.
    :param fsrc2: soft source at time n-1, (Pa).
    :type fsrc2: 2D numpy array of floats.
    :param Nb: boundary of the domain (1 if BC, 0 else).
    :type Nb: 2D numpy array of integers.
    :param c: sound speed (m.s-1).
    :type c: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param Ts: time step (s).
    :type Ts: float
    :param dx: spatial step (m).
    :type dx: float
    :param Cn: Courant number.
    :type Cn: float
    :param A: inertance of the boundary.
    :type A: float
    :param B: stiffness of the boundary.
    :type B: float
    :param C: resistivity of the boundary.
    :type C: float
    :param x_in_idx: x coordinates of the p=0 cells under the slope BC cells
    :type x_in_idx: list of integers
    :param y_in_idx: y coordinates of the p=0 cells under the slope BC cells
    :type y_in_idx: list of integers
    :param x_edges_idx: x coordinates of the slope BC cells that have single
    face contact, i.e. edges
    :type x_edges_idx: list of integers
    :param y_edges_idx: y coordinates of the slope BC cells that have single
    face contact, i.e. edges
    :type y_edges_idx: list of integers
    :param x_corners_idx: x coordinates of the slope BC cells that have two
    faces contact, i.e. corners
    :type x_corners_idx: list of integers
    :param y_corners_idx: y coordinates of the slope BC cells that have two
    faces contact, i.e. corners
    :type y_corners_idx: list of integers
    :param slope_start: index of the slope start along the x axis
    :type slope_start: int
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    j = 0
    p[:, j] = 0
    j = p.shape[1] - 1
    p[:, j] = 0
    i = 0
    p[i, :] = 0
    i = p.shape[0] - 1
    p[i, :] = 0
    # =========================================================================
    #   main solver (inside the domain): SRL FDTD <==> Cubic cell FVTD
    # =========================================================================
    A_corr_src = -1. / 3.1739566e3
    p[1:-1, 1:-1] = \
    1. / (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) + (.5 * B)))) * \
    (
        (
            (Cn ** 2) *
            (p1[2:, 1:-1] + p1[:-2, 1:-1] +
             p1[1:-1, 2:] + p1[1:-1, :-2] -
             ((4. - Nb[1:-1, 1:-1]) * p1[1:-1, 1:-1]))
        )
        + (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
        p1[1:-1, 1:-1]
        - (1. + (Nb[1:-1, 1:-1] * Cn * ((A / Ts) - (.5 * B)))) *
        p2[1:-1, 1:-1]
    ) - A_corr_src * rho * (c ** 2 * Ts ** 2 / dx ** 2) * \
        (fsrc[1:-1, 1:-1] - fsrc2[1:-1, 1:-1]) / (2 * Ts)

    # =========================================================================
    #   Zero pressure for the cells inside the cylinder
    # =========================================================================
    for ndx_0 in range(1,len(x_in_idx)):
        p[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0

    # =========================================================================
    #   Hamilton's update in DAFx 2014 Eq. (39) with "full cell" implementation
    # =========================================================================
    #   Edge-cell
    for ndx_co in range(len(x_edges_idx)):
        i = x_edges_idx[ndx_co]
        j = y_edges_idx[ndx_co]
        p[i, j] = \
        1. / (1. + (Nb[i, j] * Cn * ((A / Ts) + (.5 * B)))) * \
        (   (   (Cn ** 2) *
                (p1[i + 1, j] + p1[i - 1, j] +
                 p1[i, j + 1] + p1[i, j - 1] -
                 (3. * p1[i, j])))
            + (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
            p1[i, j]
            - (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * B)))) * p2[i, j])
    #   Corner-cell
    for ndx_co in range(len(x_corners_idx)):
        i = x_corners_idx[ndx_co]
        j = y_corners_idx[ndx_co]
        p[i, j] = \
        1. / (1. + (Nb[i, j] * Cn * ((A / Ts) + (.5 * B)))) * \
        (   (   (Cn ** 2) *
                (p1[i + 1, j] + p1[i - 1, j] +
                 p1[i, j + 1] + p1[i, j - 1] -
                 (2. * p1[i, j])   )       )
            + (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * C * Ts)))) * 2. *
            p1[i, j]
            - (1. + (Nb[i, j] * Cn * ((A / Ts) - (.5 * B)))) * p2[i, j])
    # =========================================================================
    #   B.C. on the frame Dirichlet = perfectly reflecting.
    # =========================================================================
    j = 1
    p[:slope_start - 2, j] = p[:slope_start - 2, j + 1]
    j = p.shape[1] - 2
    p[:, j] = p[:, j - 1]
    i = 1
    p[i, :] = p[i + 1, :]
    i = p.shape[0] - 2
    p[i, :] = p[i - 1, :]
    return p
