# -*- coding: utf-8 -*-
##
# \file     pe_core.py
# \title    Parabolic equation calculation core.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 14 Nov.
##

import numpy as np


def pe_pade_one_one(k, beta, dx, dy, Nx, Ny):
    """
    Initialization of the matrices for the parabolic equation, in the case
    of acoustic wave propagation (Helmholtz) in an homogeneous atmosphere
    - i.e. no temperature, no wind gradients, using the split-step Pade method.
    This function corresponds to the Pade(1, 1) expansion.
    Matrices are defined using 2D numpy arrays that contain the complex coef.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param beta: **NORMALIZED** admitance used for the bound. cond. (kg.s.m-2).
    :type k: float
    :param dx: spatial step for the x directions (m).
    :type dx: float
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the updated pressure vetcor at the location i + 1
    :rtype: 2D numpy arrays of complexes
    """
    # =========================================================================
    #   Pade(1,1) coefficients - see [blairon_phd2002, Eqs.(2.15), p.44]
    # =========================================================================
    sig = 1j * k * dx
    p1 = (1. + sig) / 4.
    q1 = (1. - sig) / 4.

    # =========================================================================
    #   Calculation of the coef. - see [blairon_phd2002, Eqs.(2.18), p.45]
    # =========================================================================
    aj = q1 * (1. / (k * dy) ** 2)
    bj = 1. + q1 * (-2. / (k * dy) ** 2)
    cj = q1 * (1. / (k * dy) ** 2)

    dj = p1 * (1. / (k * dy) ** 2)
    ej = 1. + p1 * (-2. / (k * dy) ** 2)
    fj = p1 * (1. / (k * dy) ** 2)

    # =========================================================================
    #   Filling the arrays with their respective coefficients
    # =========================================================================
    Aij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    Bij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    rng = np.arange(Nx + 1)

    Aij[rng[:-1], rng[:-1] + 1] = aj
    Aij[rng, rng] = bj
    Aij[rng[1:], rng[1:] - 1] = cj

    Bij[rng[:-1], rng[:-1] + 1] = dj
    Bij[rng, rng] = ej
    Bij[rng[1:], rng[1:] - 1] = fj

    # =========================================================================
    #   Impedance boundary condition - see [blairon_phd2002, Eqs.(2.35), p.54]
    # =========================================================================
    bj_bc = bj + 2.*1j*k*beta*dy * q1 * (1. / (k * dy) ** 2)
    cj_bc = 2. * q1 * (1. / (k * dy) ** 2)
    Aij[0, 0] = bj_bc
    Aij[0, 1] = cj_bc

    ej_bc = ej + 2.*1j*k*beta*dy * p1 * (1. / (k * dy) ** 2)
    fj_bc = 2. * p1 * (1. / (k * dy) ** 2)
    Bij[0, 0] = ej_bc
    Bij[0, 1] = fj_bc

    return Aij, Bij


def pe_pade_two_two(k, beta, dx, dy, Nx, Ny):
    """
    Initialization of the matrices for the parabolic equation, in the case
    of acoustic wave propagation (Helmholtz) in an homogeneous atmosphere
    - i.e. no temperature, no wind gradients, using the split-step Pade method.
    This module calculates the main matrices according to
    **[chevret_phd1994, Eqs.(4.66)-(4.70), p.69]**, or
    **[gauvreau_phd1999, Eqs.(3.56)-(3.59), p.98]**.
    Matrices are defined using 2D numpy arrays that contain the complex coef.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param beta: **NORMALIZED** admitance used for the bound. cond. (kg.s.m-2).
    :type k: float
    :param dx: spatial step for the x directions (m).
    :type dx: float
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the updated pressure vetcor at the location i + 1
    :rtype: 2D numpy arrays of complexes
    """
    # =========================================================================
    #   Pade coefficients RHS, LHS=conjugate, see
    #   [chevret_phd1994, Eqs.(4.63), p.68]
    # =========================================================================
    sig = 1j * k * dx
    p1 = (3. + sig) / 4.
    p2 = (sig**2 + 6.*sig + 3) / 48.

    # =========================================================================
    #   Calculation of the coefficients, see
    #   [chevret_phd1994, Eqs.(4.66)-(4.70), p.69]
    # =========================================================================
    aj = 1. + p1 * (-2. / (k * dy) ** 2) + p2 * (6. / (k * dy) ** 4)
    bj = p1 * (1. / (k * dy) ** 2) + p2 * (-4. / (k * dy) ** 2)
    cj = p2 * (1. / (k * dy) ** 4)

    aj_conj = np.conjugate(aj)
    bj_conj = np.conjugate(bj)
    cj_conj = np.conjugate(cj)

    # =========================================================================
    #   Filling the arrays with their respective coefficients
    # =========================================================================
    Aij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    Bij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    rng = np.arange(Nx + 1)

    Aij[rng[2:], rng[2:] - 2] = cj_conj
    Aij[rng[1:], rng[1:] - 1] = bj_conj
    Aij[rng, rng] = aj_conj
    Aij[rng[:-1], rng[:-1] + 1] = bj_conj
    Aij[rng[:-2], rng[:-2] + 2] = cj_conj

    Bij[rng[2:], rng[2:] - 2] = cj
    Bij[rng[1:], rng[1:] - 1] = bj
    Bij[rng, rng] = aj
    Bij[rng[:-1], rng[:-1] + 1] = bj
    Bij[rng[:-2], rng[:-2] + 2] = cj

    # =========================================================================
    #   Impedance boundary condition - see [chevret_phd1994, Eqs.(4.70), p.71]
    # =========================================================================
    aj_bc = aj + 2.*1j*k*beta*dy*bj + 4.*1j*k*beta*dy*cj
    bj_bc = bj + 2*1j*k*beta*dy*cj

    Bij[0, 0] = aj_bc
    Bij[1, 0] = bj_bc
    Bij[0, 1] = 2. * bj
    Bij[0, 2] = 2. * cj
    Bij[1, 1] = aj + cj

    aj_bc_conj = aj_conj + 2. * 1j * k * beta * dy * bj_conj + \
                 4. * 1j * k * beta * dy * cj_conj
    bj_bc_conj = bj_conj + 2 * 1j * k * beta * dy * cj_conj

    Aij[0, 0] = aj_bc_conj
    Aij[1, 0] = bj_bc_conj
    Aij[0, 1] = 2. * bj_conj
    Aij[0, 2] = 2. * cj_conj
    Aij[1, 1] = aj_conj + cj_conj

    return Aij, Bij


def pe_solver(Aij, Bij, pi):
    """
    Parabolic equation solver that calculates the matrix product and the matrix
    inversion for each spatial iteration follwing the x direction (i +1 <- i).

    :param Aij: matrix from the LHS of **[chevret_phd1994, Eqs.(4.69), p.70]**.
    :type Aij: 2darray of complexes
    :param Bij: matrix from the RHS of **[chevret_phd1994, Eqs.(4.69), p.70]**.
    :type Bij: 2darray of complexes
    :param pi: pressure at row i (Pa).
    :type pi: 1darrays of complexes (vector)
    :return: pressure at row i + 1 (Pa).
    :rtype: 1darray of complexes (vector)
    """
    # =========================================================================
    #   Calculating the pressure at row i + 1
    # =========================================================================
    return np.dot(np.linalg.inv(Aij), np.dot(Bij, pi))


def pe_scheme_pade_22_full(k, beta, dx, dy, pi, Nx, Ny):
    """
    It contains both the array filling and the solver in the same function.
    Parabolic equation scheme that solves the acoustic wave propagation eq.
    (Helmholtz) in an homogeneous atmosphere - i.e. no temperature, no wind
    gradients, using the split-step Pade method.
    This module calculates the main matrices according to
    **[chevret_phd1994, Eqs.(4.66)-(4.70), p.69]**, or
    **[gauvreau_phd1999, Eqs.(3.56)-(3.59), p.98]**.
    Matrices are defined using 2D numpy arrays that contain the complex coef.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param beta: **NORMALIZED** admitance used for the bound. cond. (kg.s.m-2).
    :type k: float
    :param dx: spatial step for the x directions (m).
    :type dx: float
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param pi: pressure at row i (Pa).
    :type pi: 1darrays of complexes (vector)
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the updated pressure vetcor at the location i + 1
    :rtype: 2D numpy arrays of complexes
    """
    # =========================================================================
    #   Pade coefficients RHS, LHS=conjugate,
    #   see [chevret_phd1994, Eqs.(4.63), p.68]
    # =========================================================================
    sig = 1j * k * dx
    p1 = (3. + sig) / 4.
    p2 = (sig**2 + 6.*sig + 3) / 48.

    # =========================================================================
    #   Calculation of the coefficients,
    #   see [chevret_phd1994, Eqs.(4.66)-(4.70), p.69]
    # =========================================================================
    aj = 1. + p1 * (-2. / (k * dy) ** 2) + p2 * (6. / (k * dy) ** 4)
    bj = p1 * (1. / (k * dy) ** 2) + p2 * (-4. / (k * dy) ** 2)
    cj = p2 * (1. / (k * dy) ** 4)

    aj_conj = np.conjugate(aj)
    bj_conj = np.conjugate(bj)
    cj_conj = np.conjugate(cj)

    # =========================================================================
    #   Filling the arrays with their respective coefficients
    # =========================================================================
    Aij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    Bij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    rng = np.arange(Nx + 1)
    Aij[rng[2:], rng[2:] - 2] = cj_conj
    Aij[rng[1:], rng[1:] - 1] = bj_conj
    Aij[rng, rng] = aj_conj
    Aij[rng[:-1], rng[:-1] + 1] = bj_conj
    Aij[rng[:-2], rng[:-2] + 2] = cj_conj
    Bij[rng[2:], rng[2:] - 2] = cj
    Bij[rng[1:], rng[1:] - 1] = bj
    Bij[rng, rng] = aj
    Bij[rng[:-1], rng[:-1] + 1] = bj
    Bij[rng[:-2], rng[:-2] + 2] = cj

    # =========================================================================
    #   Impedance boundary condition - [chevret_phd1994, Eqs.(4.70), p.71]
    # =========================================================================
    aj_bc = aj + 2.*1j*k*beta*dy*bj + 4.*1j*k*beta*dy*cj
    bj_bc = bj + 2*1j*k*beta*dy*cj

    Bij[0, 0] = aj_bc
    Bij[1, 0] = bj_bc
    Bij[0, 1] = 2. * bj
    Bij[0, 2] = 2. * cj
    Bij[1, 1] = aj + cj

    aj_bc_conj = aj_conj + 2. * 1j * k * beta * dy * bj_conj + \
                 4. * 1j * k * beta * dy * cj_conj
    bj_bc_conj = bj_conj + 2 * 1j * k * beta * dy * cj_conj

    Aij[0, 0] = aj_bc_conj
    Aij[1, 0] = bj_bc_conj
    Aij[0, 1] = 2. * bj_conj
    Aij[0, 2] = 2. * cj_conj
    Aij[1, 1] = aj_conj + cj_conj

    # =========================================================================
    #   Calculating the updated vector at row i + 1
    # =========================================================================
    return np.dot(np.linalg.inv(Aij), np.dot(Bij, pi))


def pe_scheme_pade_1_1(k, beta, dx, dy, p_ij, Nx, Ny):
    """
    It contains both the array filling and the solver in the same function.
    Parabolic equation scheme that solves the acoustic wave propagation eq.
    (Helmholtz) in an homogeneous atmosphere - i.e. no temperature, no wind
    gradients, using the split-step Pade method order (1, 1).
    This module calculates the main matrices according to
    **[blairon_phd2002, Eqs.(2.18)-(2.19), p.45]**.
    Matrices are defined using 2D numpy arrays that contain the complex coef.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param beta: **NORMALIZED** admitance used for the bound. cond. (kg.s.m-2).
    :type k: float
    :param dx: spatial step for the x directions (m).
    :type dx: float
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param p_ij: pressure at the discrete location i,j ~ (x, y) (Pa).
    :type p_ij: 2D numpy arrays of complexes
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the updated pressure array
    :rtype: 2D numpy arrays of complexes
    """
    # =========================================================================
    #   Calculation of the coefficients, see
    #   [blairon_phd2002, Eqs.(2.18)-(2.19), p.45]
    # =========================================================================
    aj = (1. - 1j * k * dx) / (4. * k ** 2 * dy ** 2)
    bj = -2. / (k * dy) ** 2 * (1. - 1j * k * dx) / 4. + 1.
    cj = aj
    dj = (1. + 1j * k * dx) / (4. * k ** 2 * dy ** 2)
    ej = -2. / (k * dy) ** 2 * (1. + 1j * k * dx) / 4. + 1.
    fj = dj

    # =========================================================================
    #   Filling the arrays with their respective coefficients
    # =========================================================================
    Aij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    Bij = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    rng = np.arange(Nx + 1)
    Aij[rng, rng - 1] = aj
    Aij[rng, rng] = bj
    Aij[rng[:-1], rng[:-1] + 1] = cj
    Bij[rng, rng - 1] = dj
    Bij[rng, rng] = ej
    Bij[rng[:-1], rng[:-1] + 1] = fj

    # =========================================================================
    #   Impedance boundary condition - [blairon_phd2002, Eqs.(2.35), p.54]
    # =========================================================================
    bj_bc = - 2. / (k * dy) ** 2 * (1. - 1j * k * dx) / 4. + 1. \
            + 2. * 1j * k * beta * dy * (1 - 1j*k*dx) / (4. * k**2 * dy**2)
    cj_bc = 2. * (1 - 1j*k*dx) / (4. * k**2 * dy**2)
    Aij[0, 0] = bj_bc
    Aij[1, 0] = cj_bc

    ej_bc = - 2. / (k * dy) ** 2 * (1. + 1j * k * dx) / 4. + 1. \
            + 2. * 1j * k * beta * dy * (1 + 1j*k*dx) / (4. * k**2 * dy**2)
    fj_bc = 2. * (1 + 1j*k*dx) / (4. * k**2 * dy**2)
    Bij[0, 0] = ej_bc
    Bij[1, 0] = fj_bc

    # =========================================================================
    #   Calculating the updated matrix using the inversion of Aij
    # =========================================================================
    return np.dot(np.linalg.inv(Aij), np.dot(Bij, p_ij))
