# -*- coding: utf-8 -*-
##
# \file     upd_fdtd.py
# \title    FDTD updates that return the pressure at each time iteration.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 6 Feb.
##

import numpy as np


def upd_p_fdtd_srl(p, p1, p2, fsrc, fsrc1, fsrc2, Nb,
                   c, rho, Ts, dx, Cn, A, B, C, depth):
    """
    This FDTD update calculates the pressure at the discrete time n+1 (p),
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
    :param depth: index of depth of the velocity and pressure layer used for
    the impedance implementation
    :type depth: int
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    j = 0;              p[:,j] = 0
    j = p.shape[1]-1;   p[:,j] = 0
    i = 0;              p[i,:] = 0
    i = p.shape[0]-1;   p[i,:] = 0
    # =========================================================================
    #   main solver (inside the domain): SRL FDTD <==> Cubic cell FVTD
    # =========================================================================
    A_corr_src = -1./3.1739566e3
    p[1:-1, depth:-1] = 1./(1.+(Nb[1:-1, depth:-1]*Cn*((A/Ts)+(.5*B)))) * \
                (
                (
                (Cn**2) *
                (p1[2:, depth:-1]+p1[:-2, depth:-1] +
                 p1[1:-1, depth + 1:]+p1[1:-1, depth - 1:-2] -
                ((4.-Nb[1:-1, depth:-1])*p1[1:-1, depth:-1]))
                )
                + (1.+(Nb[1:-1, depth:-1]*Cn*((A/Ts)-(.5*C*Ts))))*2. *
                p1[1:-1, depth:-1]
                - (1.+(Nb[1:-1, depth:-1]*Cn*((A/Ts)-(.5*B)))) *
                p2[1:-1, depth:-1]
                ) - A_corr_src * rho * (c**2*Ts**2/dx**2) *\
                (fsrc[1:-1, depth:-1] - fsrc2[1:-1, depth:-1]) / (2.*Ts)

    # =========================================================================
    #   B.C. on the frame Dirichlet = perfectly reflecting.
    # =========================================================================
    j = 1;              p[:, j] = p[:, j+1]
    j = p.shape[1]-2;   p[:, j] = p[:, j-1]
    i = 1;              p[i, :] = p[i+1, :]
    i = p.shape[0]-2;   p[i, :] = p[i-1, :]
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


def upd_vel_pbc_fdtd(   p, c, rho, Ts, dx,
                        v_x, v_y, p_bc, p1_bc,
                        K, a_k_cor, gamma_k, psi_k, psi1_k, depth):
    """
    Velocity and pressure update defined for the impedance implementation.
    The velocity normal to the boundary is required for the definition of the
    ground impedance.
    The stencil follows a Standard RectiLinear (SRL) implementation.

    :param p: updated pressure at time n+1, (Pa).
    :type p: 2D numpy array of floats.
    :param c: sound speed (m.s-1).
    :type c: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param Ts: time step (s).
    :type Ts: float
    :param dx: spatial step (m).
    :type dx: float
    :param v_x: velocity following the x-axis direction at n+1 (m.s-1)
    :type v_x: 2D numpy array of floats.
    :param v_y: velocity following the y-axis direction at n+1 (m.s-1)
    :type v_y: 2D numpy array of floats.
    :param p_bc: presure on the boundary at n+1 (Pa)
    :type p_bc: 2D numpy array of floats.
    :param p1_bc: presure on the boundary at n (Pa)
    :type p1_bc: 2D numpy array of floats.
    :param K: order of the recursive convolution method for the ground imp.
    :type K: int
    :param a_k_cor: Miki's model correction for TDBC implementation
    - see [Guillaume_jsv2011, Sec.3.2.2] and get_imped_coefts.py, line 91.
    :type a_k_cor: list of floats
    :param gamma_k: poles of the partial fraction expansion
    :type gamma_k: list of floats
    :param psi_k: impedance accumulator at time n+1 for the recursive convolut.
    :type psi_k: 2D numpy array of floats.
    :param psi1_k: impedance accumulator at time n for the recursive convolut.
    :type psi1_k: 2D numpy array of floats.
    :param depth: index of depth of the velocity and pressure layer used
    for the impedance implementation
    :type depth: int

    :param v1_x: velocity following the x-axis direction at n (m.s-1)
    :type v1_x: 2D numpy array of floats.
    :param v1_y: velocity following the y-axis direction at n (m.s-1)
    :type v1_y: 2D numpy array of floats.

    :return: the updated pressure at the time n+1 - p ; the pressure
    at the boundary - p_bc ; the velocity components -
    v_x and v_y ; and the accumulator - psi_k.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   Pressure and velocity update:
    #   enables the access to the normal velocity v_y for the BC
    # =========================================================================
    # TODO: currently unstable, see if this BC can be stable without RK-4 upd

#   Passing the pressure values to the vel-pres update (at the limit)
    p_bc[:, depth] = np.copy(p[:, depth])
#   Velocities for each direction
    v_x[1:-1, 1:depth] = v_x[1:-1, 1:depth] - \
                         Ts / (rho*dx)*(p1_bc[2:, 1:depth] -
                                        p1_bc[1:-1, 1:depth])
    v_y[1:-1, 1:depth] = v_y[1:-1, 1:depth] - \
                         Ts / (rho*dx)*(p1_bc[1:-1, 2:depth + 1] -
                                        p1_bc[1:-1, 1:depth])
#   Pressure in the bc layer
    p_bc[1:-1, 1:depth] = p1_bc[1:-1, 1:depth] - rho*c**2*Ts/dx*(
                          v_x[1:-1, 1:depth] - v_x[:-2, 1:depth] +
                          v_y[1:-1, 1:depth] - v_y[1:-1, :depth - 1])

    # =========================================================================
    #   Boundary condition (TDBC):
    #   Miki via recursive convolution (poles&coefs.: a_k, gamma_k)
    # =========================================================================
    sum_k = np.zeros((p.shape[0], p.shape[1]), dtype=np.complex128)
    j = 1   # location of the boundary on the y-axis (j coordinate)
    for k in range(K):
        exp_gamkdt = np.exp(-gamma_k[k]*Ts)
        #   Accumulator psi_k from [cotte_aiaa2009, Eq.(20)] or
        #   [guillaume_jsv2011, Eq.(18)]
        psi_k[:, j, k] = v_y[:, j]*(1 - exp_gamkdt)/gamma_k[k] + \
                         exp_gamkdt*psi1_k[:, j, k]
        # ostashev_jasa2007
        # psi_k[:, j, k] = v_y[:, j]*Ts + exp_gamkdt*psi1_k[:, j, k]
        sum_k[:, j] += a_k_cor[k]*psi1_k[:, j, k]

    #   Boundary condition from [cotte_aiaa2009, Eq.(27)] or
    #   [guillaume_jsv2011, Eq.(17)]
    p_bc[:, j] = rho * c * (v_y[:, j] + sum_k[:, j])
    # p_bc[:, j] = p_bc[:, j+1]   # reflecting
    p[:, 1:depth] = np.copy(p_bc[:, 1:depth])

    # =========================================================================
    # B.C. on the frame: Dirichlet = perfectly reflecting.
    # =========================================================================
    # Ground pressure condition may be off to avoid conflict
    # j = 1;              p_bc[:,j] = p_bc[:,j+1]
    j = p.shape[1]-2;   p_bc[:, j] = p_bc[:, j-1]
    i = 1;              p_bc[i, :] = p_bc[i+1, :]
    i = p.shape[0]-2;   p_bc[i, :] = p_bc[i-1, :]
    return p, p_bc, v_x, v_y, psi_k
