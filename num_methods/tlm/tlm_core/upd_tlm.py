# -*- coding: utf-8 -*-
##
# \file     upd_tlm.py
# \title    Updated TLM scheme that gives the pressure at each time iteration.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 14 Feb.
##
import numpy as np


def upd_p_tlm_srl_rigid(I_t, I_b, I_l, I_r):
    """
    Simplified version of the TLM calculation, good for free field propagation
    inside a rigid frame.
    Diffusion matrix, connexion laws and boundary condition for 2D TLM network.

    :param I_t: incident pulse from the top
    :type I_t: 2D numpy array of floats.
    :param I_b: incident pulse from the bottom
    :type I_b: 2D numpy array of floats.
    :param I_l: incident pulse from the left
    :type I_l: 2D numpy array of floats.
    :param I_r: incident pulse from the right
    :type I_r: 2D numpy array of floats.
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   Diffusion matrix
    # =========================================================================
    S_t = .5 * (-I_t + I_b + I_l + I_r)
    S_b = .5 * (I_t - I_b + I_l + I_r)
    S_l = .5 * (I_t + I_b - I_l + I_r)
    S_r = .5 * (I_t + I_b + I_l - I_r)

    p = 0.5 * (I_t + I_b + I_l + I_r)

    # =========================================================================
    #   Connexion laws
    # =========================================================================
    I_t[1:-1, 1:-1] = S_b[1:-1, 2:]
    I_b[1:-1, 1:-1] = S_t[1:-1, :-2]
    I_l[1:-1, 1:-1] = S_r[:-2, 1:-1]
    I_r[1:-1, 1:-1] = S_l[2:, 1:-1]

    # =========================================================================
    #   Boundary conditions (reflecting frame)
    # =========================================================================
    I_t[1:-1, -2] = S_t[1:-1, -3]
    I_b[1:-1, 1] = S_b[1:-1, 2]  # ground
    I_l[1, 1:-1] = S_l[2, 1:-1]
    I_r[-2, 1:-1] = S_r[-3, 1:-1]
    return p


def upd_p_tlm_srl(Ts, Z_air, Z_tlm, K, a_k_cor, gamma_k,
                  I_t, I_b, I_l, I_r, psi_k, psi1_k):
    """
    TLM update with the impedance boundary condition shown in
    **[guillaume_jsv2011]**.
    Diffusion matrix, connexion laws and boundary condition for 2D TLM network.

    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param Ts: time step (s).
    :type Ts: float
    :param Z_air: impedance of air
    :type Z_air: complex
    :param Z_tlm: impedance of TLM branch
    :type Z_tlm: complex
    :param Z_ratio: impedance ratio between tlm branch and air
    :type Z_ratio: complex
    :param K: order of the recursive convolution method for the ground imped.
    :type K: int
    :param a_k: residuals of the partial fraction expansion
    :type a_k: list of floats
    :param gamma_k: poles of the partial fraction expansion
    :type gamma_k: list of floats
    :param I_t: incident pulse from the top
    :type I_t: 2D numpy array of floats.
    :param I_b: incident pulse from the bottom
    :type I_b: 2D numpy array of floats.
    :param I_l: incident pulse from the left
    :type I_l: 2D numpy array of floats.
    :param I_r: incident pulse from the right
    :type I_r: 2D numpy array of floats.
    :param S_t: scattered pulse from the top
    :type S_t: 2D numpy array of floats.
    :param S_b: scattered pulse from the bottom
    :type S_b: 2D numpy array of floats.
    :param S_l: scattered pulse from the left
    :type S_l: 2D numpy array of floats.
    :param S_r: scattered pulse from the right
    :type S_r: 2D numpy array of floats.
    :param psi_k: impedance accumulator at time n+1 for the recursive convol.
    :type psi_k: 2D numpy array of floats.
    :param psi1_k: impedance accumulator at time n for the recursive convol.
    :type psi1_k: 2D numpy array of floats.
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   Diffusion matrix
    # =========================================================================
    S_t      = .5 * (-I_t + I_b + I_l + I_r)
    S_b      = .5 * (I_t - I_b + I_l + I_r)
    S_l      = .5 * (I_t + I_b - I_l + I_r)
    S_r      = .5 * (I_t + I_b + I_l - I_r)
    
    p = 0.5 * (I_t + I_b + I_l + I_r)

    # =========================================================================
    #   Connexion laws
    # =========================================================================
    I_t[1:-1, 1:-1] = S_b[1:-1, 2:]
    I_b[1:-1, 1:-1] = S_t[1:-1, :-2]
    I_l[1:-1, 1:-1] = S_r[:-2, 1:-1]
    I_r[1:-1, 1:-1] = S_l[2:, 1:-1]

    # =========================================================================
    #   Boundary conditions (reflecting frame)
    # =========================================================================
    I_t[1:-1,-2]    = S_t[1:-1,-3]
    I_b[1:-1,1]     = S_b[1:-1,1]  #   ground
    I_l[1,1:-1]     = S_l[2,1:-1]
    I_r[-2,1:-1]    = S_r[-3,1:-1]

    # =========================================================================
    #   Boundary condition (TDBC):
    #   Miki via recursive convolution (poles&coef: a_k, gamma_k)
    # =========================================================================
    # TODO: fix this TDBC in the TLM update: upd_p_tlm_srl !
    sum_k_1 = np.zeros((p.shape[0], p.shape[1]), dtype=np.complex128)
    sum_k_2 = np.zeros((p.shape[0], p.shape[1]), dtype=np.complex128)
    GAM_k = np.zeros((p.shape[0], p.shape[1]), dtype=np.complex128)

    j = 2   # location of the boundary on the y-axis (j coordinate)
    for k in range(K):
        exp_gamkdt = np.exp(-gamma_k[k]*Ts)
        #   Accumulator psi_k from [guillaume_jsv2011, Eq.(23)]
        psi_k[:, j, k] = (S_b[:, j] - S_t[:, j - 1])/Z_tlm * \
                         (1 - exp_gamkdt)/gamma_k[k] + \
                         exp_gamkdt*psi1_k[:, j, k]
        sum_k_1[:, j] += a_k_cor[k]*psi_k[:, j, k]
        sum_k_2[:, j] += a_k_cor[k]*exp_gamkdt*psi_k[:, j, k]

    # [guillaume_jsv2011, Eq.(22)]
    GAM_k[:, j] = Z_air/Z_tlm*(1. + sum_k_1[:, j])

    # [guillaume_jsv2011, Eq.(21)]
    S_t[:, j - 1] = S_b[:, j]*(-1. + GAM_k[:, j])/(1. + GAM_k[:, j]) + \
                    Z_air/(1. + GAM_k[:, j])*sum_k_2[:, j]

    I_b[1:-1, j] = S_t[1:-1, j - 1]
    return p, psi_k


def upd_p_tlm_srl_2D_scat(I_t, I_b, I_l, I_r, x_in_idx, y_in_idx,
                          x_edges_idx_1st, y_edges_idx_1st,
                          x_edges_idx_2nd, y_edges_idx_2nd,
                          x_edges_idx_3rd, y_edges_idx_3rd,
                          x_edges_idx_4th, y_edges_idx_4th,
                          x_corners_idx_1st, y_corners_idx_1st,
                          x_corners_idx_2nd, y_corners_idx_2nd,
                          x_corners_idx_3rd, y_corners_idx_3rd,
                          x_corners_idx_4th, y_corners_idx_4th):
    """
    TLM update used for case4 : scattering, in which the coordinates of the
    boundary nodes (the circular obstacle) are required.
    Diffusion matrix, connexion laws and boundary condition for 2D TLM network.

    :param I_t: incident pulse from the top
    :type I_t: 2D numpy array of floats.
    :param I_b: incident pulse from the bottom
    :type I_b: 2D numpy array of floats.
    :param I_l: incident pulse from the left
    :type I_l: 2D numpy array of floats.
    :param I_r: incident pulse from the right
    :type I_r: 2D numpy array of floats.
    :param x_in_idx: x coordinates of the p=0 cells under the slope BC cells
    :type x_in_idx: list of integers
    :param y_in_idx: y coordinates of the p=0 cells under the slope BC cells
    :type y_in_idx: list of integers
    :param x_edges_idx_1st: x coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 1st quadrant
    :type x_edges_idx_1st: list of integers
    :param y_edges_idx_1st: y coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 1st quadrant
    :type y_edges_idx_1st: list of integers
    :param x_edges_idx_2nd: x coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 2nd quadrant
    :type x_edges_idx_2nd: list of integers
    :param y_edges_idx_2nd: y coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 2nd quadrant
    :type y_edges_idx_2nd: list of integers
    :param x_edges_idx_3rd: x coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 3rd quadrant
    :type x_edges_idx_3rd: list of integers
    :param y_edges_idx_3rd: y coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 3rd quadrant
    :type y_edges_idx_3rd: list of integers
    :param x_edges_idx_4th: x coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 4th quadrant
    :type x_edges_idx_4th: list of integers
    :param y_edges_idx_4th: y coordinates of the slope BC cells that have a
    single face contact, i.e. edges, 4th quadrant
    :type y_edges_idx_4th: list of integers
    :param x_corners_idx_1st: x coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 1st quadrant
    :type x_corners_idx_1st: list of integers
    :param y_corners_idx_1st: y coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 1st quadrant
    :type y_corners_idx_1st: list of integers
    :param x_corners_idx_2nd: x coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 2nd quadrant
    :type x_corners_idx_2nd: list of integers
    :param y_corners_idx_2nd: y coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 2nd quadrant
    :type y_corners_idx_2nd: list of integers
    :param x_corners_idx_3rd: x coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 3rd quadrant
    :type x_corners_idx_3rd: list of integers
    :param y_corners_idx_3rd: y coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 3rd quadrant
    :type y_corners_idx_3rd: list of integers
    :param x_corners_idx_4th: x coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 4th quadrant
    :type x_corners_idx_4th: list of integers
    :param y_corners_idx_4th: y coordinates of the slope BC cells that have 2
    faces contact, i.e. corners, 4th quadrant
    :type y_corners_idx_4th: list of integers
    :return: the updated pressure at the time n+1, variable p.
    :rtype: 2D numpy array of floats.
    """
    # =========================================================================
    #   Diffusion matrix
    # =========================================================================
    S_t = .5 * (-I_t + I_b + I_l + I_r)
    S_b = .5 * (I_t - I_b + I_l + I_r)
    S_l = .5 * (I_t + I_b - I_l + I_r)
    S_r = .5 * (I_t + I_b + I_l - I_r)

    p = 0.5 * (I_t + I_b + I_l + I_r)

    # ========================================================================
    #   Connexion laws
    # =========================================================================
    I_t[1:-1, 1:-1] = S_b[1:-1, 2:]
    I_b[1:-1, 1:-1] = S_t[1:-1, :-2]
    I_l[1:-1, 1:-1] = S_r[:-2, 1:-1]
    I_r[1:-1, 1:-1] = S_l[2:, 1:-1]

    # =========================================================================
    #   Zero pressure for the cells inside the cylinder
    # =========================================================================
    for ndx_0 in range(len(x_in_idx)):
        I_t[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_b[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_l[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_r[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0

    # CORNERS
    for ndx_co in range(len(x_corners_idx_1st)):
        i = x_corners_idx_1st[ndx_co]
        j = y_corners_idx_1st[ndx_co]
        I_b[i, j] = S_b[i, j + 1]
        I_l[i, j] = S_l[i + 1, j]

    for ndx_co in range(len(x_corners_idx_2nd)):
        i = x_corners_idx_2nd[ndx_co]
        j = y_corners_idx_2nd[ndx_co]
        I_r[i, j] = S_r[i - 1, j]
        I_b[i, j] = S_b[i, j + 1]

    for ndx_co in range(len(x_corners_idx_3rd)):
        i = x_corners_idx_3rd[ndx_co]
        j = y_corners_idx_3rd[ndx_co]
        I_t[i, j] = S_t[i, j - 1]
        I_r[i, j] = S_r[i - 1, j]

    for ndx_co in range(len(x_corners_idx_4th)):
        i = x_corners_idx_4th[ndx_co]
        j = y_corners_idx_4th[ndx_co]
        I_t[i, j] = S_t[i, j - 1]
        I_l[i, j] = S_l[i + 1, j]

    # EDGES
    for ndx_co in range(len(x_edges_idx_1st)):
        i = x_edges_idx_1st[ndx_co]
        j = y_edges_idx_1st[ndx_co]
        I_l[i, j] = S_l[i + 1, j]

    for ndx_co in range(len(x_edges_idx_2nd)):
        i = x_edges_idx_2nd[ndx_co]
        j = y_edges_idx_2nd[ndx_co]
        I_b[i, j] = S_b[i, j + 1]

    for ndx_co in range(len(x_edges_idx_3rd)):
        i = x_edges_idx_3rd[ndx_co]
        j = y_edges_idx_3rd[ndx_co]
        I_r[i, j] = S_r[i - 1, j]

    for ndx_co in range(len(x_edges_idx_4th)):
        i = x_edges_idx_4th[ndx_co]
        j = y_edges_idx_4th[ndx_co]
        I_t[i, j] = S_t[i, j - 1]

    # =========================================================================
    #   Boundary conditions (reflecting frame)
    # =========================================================================
    I_t[1:-1, -2] = S_t[1:-1, -3]
    I_b[1:-1, 1] = S_b[1:-1, 2]  # ground
    I_l[1, 1:-1] = S_l[2, 1:-1]
    I_r[-2, 1:-1] = S_r[-3, 1:-1]
    return p

def upd_p_tlm_srl_2D_slope(I_t, I_b, I_l, I_r,
                           x_in_idx, y_in_idx, x_edges_idx, y_edges_idx,
                           x_corners_idx, y_corners_idx, slope_start_idx):

    """
    TLM update for case5: horizontal plane ended by an upward sloping part.
    Diffusion matrix, connexion laws and boundary condition for 2D TLM network.

    :param I_t: incident pulse from the top
    :type I_t: 2D numpy array of floats.
    :param I_b: incident pulse from the bottom
    :type I_b: 2D numpy array of floats.
    :param I_l: incident pulse from the left
    :type I_l: 2D numpy array of floats.
    :param I_r: incident pulse from the right
    :type I_r: 2D numpy array of floats.
    :param x_in_idx: x coordinates of the p=0 cells under the slope BC cells
    :type x_in_idx: list of integers
    :param y_in_idx: y coordinates of the p=0 cells under the slope BC cells
    :type y_in_idx: list of integers
    :param x_edges_idx: x coordinates of the slope BC cells that have a
    single face contact, i.e. edges
    :type x_edges_idx: list of integers
    :param y_edges_idx: y coordinates of the slope BC cells that have a
    single face contact, i.e. edges
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
    # =========================================================================
    #   Diffusion matrix
    # =========================================================================
    S_t = .5 * (-I_t + I_b + I_l + I_r)
    S_b = .5 * (I_t - I_b + I_l + I_r)
    S_l = .5 * (I_t + I_b - I_l + I_r)
    S_r = .5 * (I_t + I_b + I_l - I_r)

    p = 0.5 * (I_t + I_b + I_l + I_r)

    # =========================================================================
    #   Connexion laws
    # =========================================================================
    I_t[1:-1, 1:-1] = S_b[1:-1, 2:]
    I_b[1:-1, 1:-1] = S_t[1:-1, :-2]
    I_l[1:-1, 1:-1] = S_r[:-2, 1:-1]
    I_r[1:-1, 1:-1] = S_l[2:, 1:-1]

    # =========================================================================
    #   Zero pressure for the cells inside the cylinder
    # =========================================================================
    for ndx_0 in range(len(x_in_idx)):
        I_t[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_b[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_l[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0
        I_r[x_in_idx[ndx_0], y_in_idx[ndx_0]] = 0

    # CORNER
    for ndx_co in range(len(x_corners_idx)):
        i = x_corners_idx[ndx_co]
        j = y_corners_idx[ndx_co]
        I_r[i, j] = S_r[i - 1, j]
        I_b[i, j] = S_b[i, j + 1]

    # EDGE
    for ndx_co in range(len(x_edges_idx)):
        i = x_edges_idx[ndx_co]
        j = y_edges_idx[ndx_co]
        I_b[i, j] = S_b[i, j + 1]

    # =========================================================================
    #   Boundary conditions (reflecting frame)
    # =========================================================================
    I_t[1:-1, -2] = S_t[1:-1, -3]
    I_b[1:slope_start_idx - 2, 1] = S_b[1:slope_start_idx - 2, 2]  # ground
    I_l[1, 1:-1] = S_l[2, 1:-1]
    I_r[-2, 1:-1] = S_r[-3, 1:-1]
    return p