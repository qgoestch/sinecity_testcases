# -*- coding: utf-8 -*-
##
# \file     init_tlm_upslope.py
# \title    Definition of the numerical parameters for the TLM method for the
#           grid convergence sudy. The satial is an input parameter.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_tlm.py files.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 06 Oct.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

tlm_core_path = os.path.join(base_path, 'tlm_core')
site.addsitedir(tlm_core_path)
from upd_tlm import upd_p_tlm_srl_2D_slope

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
import source_signals as src

#data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
#site.addsitedir(data_plotting_path)
#from plot_boundary_cells import plot_cells_slope
#from display_wavefronts import instatenous_pressure


def tlm_srl_init_slope(freq, h_s, l1max, l2max, zmax, alpha, c,
                        dl, h_num, h_set, T, rho, case, disp_inst_p):
    """
         Setting the 2D geometries and running the TLM update for case 5: slope.
      Main script that contains all the parameters to run the TLM update in 2D.

    :param freq: frequency of the source signal (Hz).
    :type freq: float
    :param h_s: source heights (m).
    :type h_s: float
    :param l1max: length of the horizontal part (m).
    :type l1max: float
    :param l2max: length of the upward-sloping part (m).
    :type l2max: float
    :param zmax: height of the the two parts (m).
    :type zmax: float
    :param alpha: angle of the upward-sloping part (rad).
    :type alpha: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param T: simulation duration (s).
    :type T: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param disp_inst_p: display the instantaneous pressure.
    :type disp_inst_p: bool

    :param      src_typ source type, string "gauss_1", "sine", "ricker" "dirac"...
    :param      src_frq source freauency, scalar (Hz).

    :param      T       simulation duration, scalar (s).
    :param      dl_max  maximal spatial step regarding the source frequency, scalar (m).
    :param      dt      time step, scalar (s).
    :param      Lx      continuous length of the domain (in meter) following the x-direction.
    :param      Ly      continuous length of the domain (in meter) following the y-direction.
    :param      Nx      discrete length of the domain (number of node) following the x-direction.
    :param      Ny      discrete length of the domain (number of node) following the y-direction.
    :param      x       discrete length sequence of a domain side, scalar (m).
    :param      dx      spatial step after discretization, scalar (m).
    :param      Nt      number of iteration time, scalar.
    :param      t       discretized time sequence, 1d array (s).
    :param      n       discrete iteration inside the for loop, scalar.
    :param      It      discretized time sequence, 1d array.
    :param      Ts      time step after dicretization, scalar (s).

    :param      src_dly source delay in number of iteration, scalar.

    :param      x_src   discrete x coordinate of the source, scalar (number of node).
    :param      y_src   discrete y coordinate of the source, scalar (number of node).
    :param      x_rcv   discrete x coordinate of the receiver, scalar (number of node).
    :param      y_rcv   discrete y coordinate of the receiver, scalar (number of node).

    :param      Z_tlm   impedance of TLM branch, scalar.
    :param      Z_air   impedance of air, scalar.
    :param      Z_ratio impedance ratio between tlm branch and air, scalar.
    :param      K       order of the recursive convolution method for the ground impedance, scalar.
    :param      a_k     residuals of the partial fraction expansion, list of K size.
    :param      gamma_k poles of the partial fraction expansion, list of K size.

    :param      I_t     incident pulse from the top.
    :param      I_b     incident pulse from the bottom.
    :param      I_l     incident pulse from the left.
    :param      I_r     incident pulse from the right.
    :param      S_t     scattered pulse from the top.
    :param      S_b     scattered pulse from the bottom.
    :param      S_l     scattered pulse from the left.
    :param      S_r     scattered pulse from the right.

    :param      psi_k   impedance accumulator at n+1 for the recursive convolution, numpy array (dimension of the scene).
    :param      psi1_k  impedance accumulator at n for the recursive convolution, numpy array (dimension of the scene).

    :param      geo     boundary of the domain (1 if BC, 0 else), similar to Nb
                        in FDTD scripts, numpy array (dimension of the scene).

    :return: the acoustic pressure at the pre-defined receivers' locations as a function of time.
    :rtype: (2+1)D array of floats
    """
    # ==============================================================================
    #   Source
    # ==============================================================================
    src_typ = "gauss_1"
    src_frq = freq

    # ==============================================================================
    #   Parameters
    # ==============================================================================
    # dl_max = 0.075 * c * np.sqrt(2) / src_frq  # 2% accuracy kowalczyk_ieee2011 ~ lambda/9.43
    dt = dl / (np.sqrt(2.) * c)
    Ly = 3. * zmax
    Lx = 4. * (l1max + l2max) + (l1max + np.cos(alpha)*l2max)
    Nx = np.int(np.round(Lx / dl))
    Ny = np.int(np.round(Ly / dl))
    x = np.linspace(0, Lx, Nx + 1)
    dx = np.float64(x[1] - x[0]);
    dx = round(dx, 5)
    Nt = int(round(T / float(dt)))
    t = np.linspace(0, Nt * dt, Nt + 1)
    It = range(0, t.shape[0])
    Ts = np.float64(t[1] - t[0])

    # Cn_lim = np.sqrt(2)**-1
    # dt_coarse = 2 * 10 ** -4
    # dl_coarse = h_set[-1]
    # c = np.float64(Cn_lim * dl_coarse / (dt_coarse))

    src_dly = int(T / 20. / Ts)

    print '                     TLM above upward slope                    '
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
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' % (Nx, Ny, Lx, Ly)
    # print 'TL-IMPEDANCE:    Z_tlm=%f ., Z_tlm_2=%f .' % (Z_tlm, Z_tlm_2)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % (src_typ)
    print 'SOURCE FREQ:  f=%g Hz.' % (src_frq)
    print '---------------------------------------------------------------'
    # ==============================================================================
    #   Incident (I) and Reflected (R) pulses from all directions (top, bottom...)
    # ==============================================================================
    I_t = np.zeros((Nx + 1, Ny + 1))
    I_b = np.zeros((Nx + 1, Ny + 1))
    I_l = np.zeros((Nx + 1, Ny + 1))
    I_r = np.zeros((Nx + 1, Ny + 1))
    S_t = np.zeros((Nx + 1, Ny + 1))
    S_b = np.zeros((Nx + 1, Ny + 1))
    S_l = np.zeros((Nx + 1, Ny + 1))
    S_r = np.zeros((Nx + 1, Ny + 1))
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    bc_ed = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
    bc_co = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
    bc_in = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)

    # ==============================================================================
    #   Boundaries of the domain, where the BC is calculated
    # ==============================================================================
    geo = np.zeros((Nx + 1, Ny + 1))
    i = 1;
    geo[i, 1:-1] = 1.
    i = p.shape[0] - 2;
    geo[i, 1:-1] = 1.
    j = 1;
    geo[1:-1, j] = 1.
    j = p.shape[1] - 2;
    geo[1:-1, j] = 1.

    # ==============================================================================
    #   Definition of the slope
    # ==============================================================================
    x_slope = np.arange(2. * (l1max + l2max) + l1max, Lx, dx)
    y_slope = np.arange(1. * dx, (Lx - (2. * (l1max + l2max) + l2max)) + 1. * dx, dx) * np.tan(alpha*np.pi/180)
    x_slope_idx = [int(round(i / dx)) for i in x_slope]
    y_slope_idx = [int(round(i / dx)) for i in y_slope]

    x_in_idx = x_slope_idx
    y_in_idx = [i - 1 for i in y_slope_idx]
    for ndx in range(len(x_in_idx)):
        bc_in[x_in_idx[ndx], y_in_idx[ndx]] = 4

    x_corner_idx = []
    y_corner_idx = []
    x_edge_idx = []
    y_edge_idx = []
    for ndx_co in range(len(x_slope_idx) - 1):
        if  ((x_slope_idx[ndx_co] == x_in_idx[ndx_co]) and
            (y_slope_idx[ndx_co] - 1 == y_in_idx[ndx_co])) and \
            ((x_slope_idx[ndx_co] + 1 == x_in_idx[ndx_co + 1]) and
            (y_slope_idx[ndx_co]  == y_in_idx[ndx_co + 1])):
            x_corner_idx.append(x_slope_idx[ndx_co])
            y_corner_idx.append(y_slope_idx[ndx_co])

        if  ((x_slope_idx[ndx_co] == x_in_idx[ndx_co]) and
            (y_slope_idx[ndx_co] - 1 == y_in_idx[ndx_co])) and \
            ((x_slope_idx[ndx_co] + 1 != x_in_idx[ndx_co + 1]) or
            (y_slope_idx[ndx_co]  != y_in_idx[ndx_co + 1])):
            x_edge_idx.append(x_slope_idx[ndx_co])
            y_edge_idx.append(y_slope_idx[ndx_co])

    # for ndx in range(len(x_edge_idx)):
    #     bc_ed[x_edge_idx[ndx], y_edge_idx[ndx]] = 1
    # for ndx in range(len(x_corner_idx)):
    #     bc_co[x_corner_idx[ndx], y_corner_idx[ndx]] = 2

    slope_start = int(round((2. * (l1max + l2max) + l1max) / dx))
    Nz_max = int(round(zmax/dx))
    slope_start_idx = int(round((2. * (l1max + l2max) + l1max) / dx))
    # plot_cells_slope(bc_ed, bc_co, bc_in, slope_start, Nx, Nz_max)

    # ==============================================================================
    #   Location of the source and receiver(s)
    # ==============================================================================
    x_src_idx = int(round((2. * (l1max + l2max)) / dx))
    y_src_idx = int(round((h_s / dx)))

    x_st_l1_idx = int(round((2. * (l1max + l2max)) / dx))
    x_en_l1_idx = x_st_l1_idx + int(round((l1max / dx)))
    y_st_l1_idx = 1
    y_en_l1_idx = int(round(zmax / dx))

    x_st_l2_idx = x_en_l1_idx + 1
    x_en_l2_idx = x_st_l2_idx + int(round((l2max / dx)))
    x_l2_seq_idx = np.arange(x_st_l2_idx,x_en_l2_idx, 1, dtype=int)
    y_l2_seq = np.arange(1. * dx, (Lx - (2. * (l1max + l2max) + l2max)) + 1. * dx, dx) * np.tan(alpha*np.pi/180)
    y_l2_seq_idx = [int(round(i / dx)) for i in y_l2_seq]

    p_saved_l1 = np.zeros((int(round(l1max / dx)), int(round(zmax / dx)) - 1, Nt), dtype=np.complex128)
    # p_saved_l2 = np.zeros((int(round(l2max / dx)), len(y_l2_seq_idx), Nt), dtype=np.complex128)
    p_saved_l2 = np.zeros((int(round(l2max / dx)), int(round(zmax / dx)) - 1, Nt), dtype=np.complex128)

    # ==============================================================================
    #   Hard source assignment
    # ==============================================================================
    A_fdtd_inv = 1.e0  # see gauss_1 in source_signals.py
    A_tlm = 5. * A_fdtd_inv
    for n in It[:-1]:
        # source signal impl.
        I_t[x_src_idx, y_src_idx] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_b[x_src_idx, y_src_idx] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_l[x_src_idx, y_src_idx] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_r[x_src_idx, y_src_idx] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)

        # ==============================================================================
        #   Calculation of the Incident and Scattered pulses (via upd_tlm.py)
        # ==============================================================================
        p = upd_p_tlm_srl_2D_slope(I_t, I_b, I_l, I_r,
                                   x_in_idx, y_in_idx, x_edge_idx, y_edge_idx,
                                   x_corner_idx, y_corner_idx, slope_start_idx)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, dt, Lx, Ly, case, True)

        p_saved_l1[:, :, n] = p[x_st_l1_idx: x_en_l1_idx, y_st_l1_idx: y_en_l1_idx]

        p_saved_l2[:, :, n] = p[x_st_l2_idx: x_en_l2_idx, y_st_l1_idx: y_en_l1_idx]

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'tlm')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % h_num), t)
    np.save(os.path.join(res_path, 'Ts_%i.npy' % h_num), Ts)
    np.savez_compressed(os.path.join(res_path, 'p_%i_%ideg.npz' % (h_num, int(alpha))),
                        p_l1=p_saved_l1, p_l2=p_saved_l2)
    # np.savez_compressed(os.path.join(res_path, 'p_l1_%i_%ideg.npz' % (h_num, int(alpha))), p_saved_l1)
    # np.savez_compressed(os.path.join(res_path, 'p_l2_%i_%ideg.npz' % (h_num, int(alpha))), p_saved_l2)
