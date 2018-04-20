# -*- coding: utf-8 -*-
##
# \file     init_tlm_ground.py
# \title    Definition of the numerical parameters for the TLM method for the
#           grid convergence sudy. The satial is an input parameter.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_tlm.py files.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 09 Aug.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

tlm_core_path = os.path.join(base_path, 'tlm_core')
site.addsitedir(tlm_core_path)
from upd_tlm import upd_p_tlm_srl_rigid, upd_p_tlm_srl

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from get_imped_coefts import get_coefts_Miki
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure


def tlm_srl_init_impgr(dt, dl, h_num, h_set, d_sr, h_s, h_r, T, f_max, rho,
                       sigma, case, free_field, disp_inst_p):
    """
    Setting the 2D geometries and running the TLM update for case 3:
    ground reflection.
    Main script that contains all the parameters to run the TLM update in 2D.
    
    :param dt: time step (s).
    :type dt: float
    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param T: simulation duration (s).
    :type T: float
    :param f_max: approximated maximale frequency of the source signal (Hz).
    :type f_max: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param sigma: pecific airflow resistivity (kNm-4s==CGS).
    :type sigma: float
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param free_field: the domain is enlarged
    :type free_field: bool
    :param disp_inst_p: display the instantaneous pressure.
    :type disp_inst_p: bool

    :param Nx: discrete length of the domain (number of node) following
    the x-axis direction.
    :param Ny: discrete length of the domain (number of node) following
    the y-axis direction.
    :param x: discrete length sequence of a domain side, scalar (m).
    :param dx: spatial step after discretization, scalar (m).
    :param Nt: number of iteration time, scalar.
    :param t: discretized time sequence, 1d array (s).
    :param n: discrete iteration inside the for loop, scalar.
    :param It: discretized time sequence, 1d array.
    :param Ts: time step after dicretization, scalar (s).

    :param src_dly: source delay in number of iteration, scalar.

    :param x_src: discrete x coordinate of the source, int (number of node).
    :param y_src: discrete y coordinate of the source, int (number of node).
    :param x_rcv: discrete x coordinate of the receiver, int (number of node).
    :param y_rcv: discrete y coordinate of the receiver, int (number of node).

    :param Z_tlm: impedance of TLM branch, scalar.
    :param Z_air: impedance of air, scalar.
    :param Z_ratio: impedance ratio between tlm branch and air, scalar.
    :param K: order of the recursive convolution method for the ground imp, int.
    :param a_k: residuals of the partial fraction expansion, list of K size.
    :param gamma_k: poles of the partial fraction expansion, list of K size.

    :param I_t: incident pulse from the top.
    :param I_b: incident pulse from the bottom.
    :param I_l: incident pulse from the left.
    :param I_r: incident pulse from the right.
    :param S_t: scattered pulse from the top.
    :param S_b: scattered pulse from the bottom.
    :param S_l: scattered pulse from the left.
    :param S_r: scattered pulse from the right.

    :param psi_k: impedance accumulator at n+1 for the recursive convolution,
    numpy array (dimension of the scene).
    :param psi1_k: impedance accumulator at n for the recursive convolution,
    numpy array (dimension of the scene).

    :param geo: boundary of the domain (1 if BC, 0 else), similar to Nb
                        in FDTD scripts, numpy array (dimension of the scene).

    :return: the acoustic pressure at the pre-defined receivers' locations as a
    function of time.
    :rtype: (1+1)D array of floats 64
    """
    # =========================================================================
    #   Source
    # =========================================================================
    src_typ = "gauss_1"
    src_frq = f_max

    # =========================================================================
    #   Parameters
    # =========================================================================
    # dl_max = 0.075 * c * np.sqrt(2) / src_frq
    #  2% accuracy kowalczyk_ieee2011 ~ lambda/9.43
    c = 340.00
    dt = dl / (np.sqrt(2.) * c)
    Lx = 3. * d_sr[-1]
    if free_field:
        Ly = 5. * max(h_s, h_r[-1])
    else:
        Ly = 2.2 * max(h_s, h_r[-1])
    Nx = np.int(np.round(Lx / dl))
    Ny = np.int(np.round(Ly / dl))
    x = np.linspace(0, Lx, Nx + 1)
    dx = np.float64(x[1] - x[0])
    dx = round(dx, 5)
    Nt = int(round(T / float(dt)))
    t = np.linspace(0, Nt * dt, Nt + 1)
    It = range(0, t.shape[0])
    Ts = np.float64(t[1] - t[0])

    Cn_lim = np.sqrt(2)**-1
    c = np.float64(Cn_lim * dx / Ts)

    src_dly = int(T / 2. / Ts)

    # =========================================================================
    #   Set the material parameters for the BC calculation (in upd_tlm)
    # =========================================================================
    c_tlm = dx / Ts
    Z_tlm = rho * c_tlm
    Z_air = rho * c
    K = 6
    a_k_cor, gamma_k, a_k = get_coefts_Miki(K, sigma)

    # =========================================================================
    #   Receiver(s) location
    # =========================================================================
    x_src = int(round(Lx / 3 / dx))
    x_rcv = [int(round(i / dx) + x_src) for i in d_sr]
    if free_field:
        # BC length correction: 2*dx above, BC at index 0
        y_src = int(round(h_s / dx) + round(Ny / 2.)) + 1
        y_rcv = [int(round(i / dx) + round(Ny / 2.)) + 1 for i in h_r]
    else:
        y_src = int(round(h_s / dx)) + 1
        y_rcv = [int(round(i / dx)) + 1 for i in h_r]

    if free_field:
        print '                     TLM in free field                         '
    else:
        print '                     TLM above a ground                        '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2e s' % (Ts)
    print 'NUMBER OF It: Nt=%i' % (Nt)
    print 'DURATION:     T=%.3e s,' % (T)
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    # print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    # print '2PERCENT ACCURACY:  %.3f Hz,' % (0.075 / Ts)
    # print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    # print 'SPATIAL-STEP: dx=%g m, dl_max=%g m.' % (dx, dl_max)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' \
          % (Nx, Ny, Lx, Ly)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % (src_typ)
    print 'SOURCE FREQ:  f=%g Hz.' % (src_frq)
    print '---------------------------------------------------------------'
    # =========================================================================
    #   Incident (I) and Reflected (R) pulses from all directions
    # =========================================================================
    I_t = np.zeros((Nx + 1, Ny + 1))
    I_b = np.zeros((Nx + 1, Ny + 1))
    I_l = np.zeros((Nx + 1, Ny + 1))
    I_r = np.zeros((Nx + 1, Ny + 1))
    S_t = np.zeros((Nx + 1, Ny + 1))
    S_b = np.zeros((Nx + 1, Ny + 1))
    S_l = np.zeros((Nx + 1, Ny + 1))
    S_r = np.zeros((Nx + 1, Ny + 1))
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p_saved = np.zeros((len(d_sr), len(h_r), Nt), dtype=np.complex128)
    psi_k = np.zeros((Nx + 1, Ny + 1, K), dtype=np.complex128)
    psi1_k = np.zeros((Nx + 1, Ny + 1, K), dtype=np.complex128)

    # =========================================================================
    #   Boundaries of the domain, where the BC is calculated
    # =========================================================================
    geo = np.zeros((Nx + 1, Ny + 1))
    i = 1;
    geo[i, 1:-1] = 1.
    i = p.shape[0] - 2;
    geo[i, 1:-1] = 1.
    j = 1;
    geo[1:-1, j] = 1.
    j = p.shape[1] - 2;
    geo[1:-1, j] = 1.

    # =========================================================================
    #   Hard source assignment
    # =========================================================================
    A_fdtd_inv = 1.e0  # see gauss_1 in source_signals.py
    A_tlm = 0.5 * A_fdtd_inv
    for n in It[:-1]:
        # source signal impl.
        I_t[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_b[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_l[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_r[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)

        # =====================================================================
        #   Calculation of the Incident and Scattered pulses (via upd_tlm.py)
        # =====================================================================
        # without impedance
        p = upd_p_tlm_srl_rigid(I_t, I_b, I_l, I_r)
        # # with impedance
        # p, psi_k = upd_p_tlm_srl(Ts, Z_air, Z_tlm, K, a_k_cor, gamma_k,
        #                          I_t, I_b, I_l, I_r, psi_k, psi1_k)
        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, False)

        for d in range(len(d_sr)):
            for h in range(len(h_r)):
                p_saved[d, h, n] = p[x_rcv[d], y_rcv[h]]

        psi1_k[:, :, :] = psi_k.copy()

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i'
                            % case, 'tlm')
    if not os.path.exists(res_path): os.makedirs(res_path)
    if free_field:
        np.save(os.path.join(res_path, 't_%i.npy' % h_num), t)
        np.save(os.path.join(res_path, 'Ts_%i.npy' % h_num), Ts)
        np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('f', h_num)), p_saved)
    else:
        np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('t', h_num)), p_saved)
