# -*- coding: utf-8 -*-
##
# \file     init_fdtd_upslope.py
# \title    Definition of the numerical parameters for the FDTD method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_fdtd.py files.
#           It is applied for the grid convergence studies diffraction by an upward sloping part.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 26 Sep.
##

import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fdtd_core')
site.addsitedir(fdtd_core_path)
from upd_fdtd import upd_p_fdtd_srl_2D_slope

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from plot_boundary_cells import plot_cells_slope
from display_wavefronts import instatenous_pressure


def fdtd_srl_init_slope(freq, h_s, l1max, l2max, zmax, alpha, c,
                        dl, h_num, h_set, T, rho, case, disp_inst_p):
    """
         Setting the 2D geometries and running the FDTD update for case 5: slope.
      Main script that contains all the parameters to run the FDTD update in 2D.


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

    :param      dl      spatial step, scalar (m).
    :param      dt      time step, scalar (s).
    :param      Lx      continuous length of the domain (in meter) following the x-direction.
    :param      Ly      continuous length of the domain (in meter) following the y-direction.
    :param      Nx      discrete length of the domain (number of node) following the x-direction.
    :param      Ny      discrete length of the domain (number of node) following the y-direction.
    :param      x       discrete length sequence of a domain side, scalar (m).
    :param      dx      spatial step after discretization, scalar (m).
    :param      Nt      number of iteration time, scalar.
    :param      t       discretized time sequence, 1d array (s).
    :param      It      discretized time sequence, 1d array.
    :param      Ts      time step after dicretization, scalar (s).
    :param      Cn      Courant number, scalar.

    :param      p       updated pressure (n+1), numpy array (dimension of the scene).
    :param      p1      current pressure (n), numpy array (dimension of the scene).
    :param      p2      past pressure (n-1), numpy array (dimension of the scene).
    :param      fsrc    soft source (n+1), numpy array (dimension of the scene).
    :param      fsrc1   soft source (n), numpy array (dimension of the scene).
    :param      p_saved pressure saved at the receiver location, 1d array (length of the time sequence).

    :param      A       inertance of the boundary, scalar.
    :param      B       stiffness of the boundary, scalar.
    :param      C       resistivity of the boundary, scalar.
    :param      Nb      boundary of the domain (1 if BC, 0 else) for the compact
                        pressure update, numpy array (dimension of the scene).

    :param      x_slope_idx     x coordinates of the slope, list of integers.
    :param      y_slope_idx     y coordinates of the slope, list of integers.
    :param      x_in_idx        x coordinates of the zero cells below the slope, list of integers.
    :param      y_in_idx        y coordinates of the zero cells below the slope, list of integers.
    :param      x_corner_idx    x coordinates of the corner - 2 bc connexion - cells of the slope, list of integers.
    :param      y_corner_idx    y coordinates of the corner - 2 bc connexion - cells of the slope, list of integers.
    :param      x_edge_idx      x coordinates of the edge cells - 1 bc connexion -  of the slope, list of integers.
    :param      y_edge_idx      y coordinates of the edge cells - 1 bc connexion - of the slope, list of integers.

    :param      slope_start     grid index at which the slope starts following the x-direction, integer.

    :param      x_src   discrete x coordinate of the source, scalar (number of node).
    :param      y_src   discrete y coordinate of the source, scalar (number of node).
    :param      x_rcv   discrete x coordinate of the receiver, scalar (number of node).
    :param      y_rcv   discrete y coordinate of the receiver, scalar (number of node).

    :param      n       discrete iteration inside the for loop, scalar.

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
    dl_max = 0.075 * c * np.sqrt(2) / src_frq  # 2% accuracy kowalczyk_ieee2011 ~ lambda/9.43
    dt = dl / (np.sqrt(2.) * c)
    Ly = 3. * zmax
    Lx = 4. * (l1max + l2max) + (l1max + np.cos(alpha)*l2max)
    Nx = np.int(np.round(Lx / dl))
    Ny = np.int(np.round(Ly / dl))
    x = np.linspace(0, Lx, Nx + 1)
    dx = np.float64(x[1] - x[0])
    dx = round(dx, 5)
    Nt = int(round(T / float(dt)))
    t = np.linspace(0, Nt * dt, Nt + 1)
    It = range(0, t.shape[0])
    Ts = np.float64(t[1] - t[0])
    # Cn_lim = np.sqrt(2)**-1
    Cn = np.float64(c * Ts / dx)

    nNt = 1
    print '--------------------- Courant Number --------------------------'
    while Cn >= (1 / np.sqrt(2)):
        print 'COURANT NUMBER CORRECTION!!!'
        nNt = 1 + nNt
        t = np.linspace(0, Nt * dt, Nt + (nNt))  # time discretization
        It = range(0, t.shape[0] - 1)  # time iterations range
        Ts = np.float64(t[1] - t[0])  # sampling period for staggered grid
        Cn = np.float64(c * Ts / dx)
    print 'Additional iterations for stab: %i' % (nNt)
    Nt = Nt + nNt  # add the additional time steps to the main time sequence
    print 'Ratio Cn/Cn_th=%g < 1.0' % (Cn * np.sqrt(2))

    src_dly = int(T / 20. / Ts)

    print '    FDTD above a ground: flat part + upward sloping part       '
    print '---------------- Geometrical config ---------------------------'
    print 'flat: l1=%.2f m ; slope: l2=%.2f m ; angle=%.2f deg.' % (l1max, l2max, alpha)
    print 'domains height=%.2f m ; source height = %.2f m' % (zmax, h_s)
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2e s' % (Ts)
    print 'NUMBER OF It: Nt=%i' % (Nt)
    print 'DURATION:     T=%.3e s,' % (T)
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    print 'Sound speed:  c =%.2f m.s-1' % c
    # print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    print '2PERCENT ACCURACY:  %.3f Hz,' % (0.075 / Ts)
    # print 'COURANT NUMBER: Cn = %.3f' % (Cn)
    # print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    # print 'SPATIAL-STEP: dx=%g m, dl_max=%g m.' % (dx, dl_max)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' % (Nx, Ny, Lx, Ly)
    #    print 'IMPEDANCE:    Yb=%f .' %(Yb)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % src_typ
    print 'SOURCE FREQ:  f=%g Hz.' % src_frq
    print 'SOURCE DELAY: %0.2e s, %i n' %(src_dly*Ts, src_dly)
    print '---------------------------------------------------------------'
    print 'dl_max=%0.3f' % dl_max
    # ==============================================================================
    #   Variables
    # ==============================================================================
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    p1 = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    p2 = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    fsrc = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    fsrc1 = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    fsrc2 = np.zeros((Nx + 1, Ny + 1), dtype=np.float32)
    bc_ed = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
    bc_co = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
    bc_in = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)

    # ==============================================================================
    #   Boundaries of the domain, where the BC is calculated
    # ==============================================================================
    A = 0.
    B = 1.0
    C = 0.
    Nb = np.zeros((Nx + 1, Ny + 1))
    i = 1
    Nb[i, 1:-1] = 1.
    i = p.shape[0] - 2
    Nb[i, 1:-1] = 1.
    j = Ny - 1
    Nb[1:-1, j] = 1.
    j = p.shape[1] - 2
    Nb[1:-1, j] = 1.

    # ==============================================================================
    #   Definition of the slope
    # ==============================================================================
    x_slope = np.arange(2. * (l1max + l2max) + l1max, Lx, dx)
    y_slope = np.arange(1. * dx, (Lx - (2. * (l1max + l2max) + l2max)) + 1. * dx, dx) * np.tan(alpha*np.pi/180)
    x_slope_idx = [int(round(i / dx)) for i in x_slope]
    y_slope_idx = [int(round(i / dx)) for i in y_slope]

    x_in_idx = x_slope_idx                      #+ x_slope_idx
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
            (y_slope_idx[ndx_co] != y_in_idx[ndx_co + 1])):
            x_edge_idx.append(x_slope_idx[ndx_co])
            y_edge_idx.append(y_slope_idx[ndx_co])

    # for ndx in range(len(x_edge_idx)):
    #     bc_ed[x_edge_idx[ndx], y_edge_idx[ndx]] = 1
    #     Nb[x_edge_idx[ndx], y_edge_idx[ndx]] = 1.
    #
    # for ndx in range(len(x_corner_idx)):
    #     bc_co[x_corner_idx[ndx], y_corner_idx[ndx]] = 2
    #     Nb[x_corner_idx[ndx], y_corner_idx[ndx]] = 2.

    # # Additional layers of zero pressure
    # x_in_idx = x_slope_idx[1:]  # + x_slope_idx
    # y_in_idx = [i - 1 for i in y_slope_idx[:-1]]
    # for ndx in range(len(x_in_idx)):
    #     bc_in[x_in_idx[ndx], y_in_idx[ndx]] = 4

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
    x_l2_seq_idx = np.arange(x_st_l2_idx, x_en_l2_idx, 1, dtype=int)
    y_l2_seq = np.arange(1. * dx, (Lx - (2. * (l1max + l2max) + l2max)) + 1. * dx, dx) * np.tan(alpha * np.pi / 180)
    y_l2_seq_idx = [int(round(i / dx)) for i in y_l2_seq]

    p_saved_l1 = np.zeros((int(round(l1max / dx)), int(round(zmax / dx)) - 1, Nt), dtype=np.complex128)
    # p_saved_l2 = np.zeros((int(round(l2max / dx)), len(y_l2_seq_idx), Nt), dtype=np.complex128)
    p_saved_l2 = np.zeros((int(round(l2max / dx)), int(round(zmax / dx)) - 1, Nt), dtype=np.complex128)

    # ==============================================================================
    #   Calculation of the pressure
    # ==============================================================================
    depth = 1
    for n in It[:]:
        # p1[x_src_idx, y_src_idx] = 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # hard source impl.
        fsrc[x_src_idx, y_src_idx] = 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # soft source impl.
        p = upd_p_fdtd_srl_2D_slope(p, p1, p2, fsrc, fsrc2, Nb, c, rho, Ts, dx, Cn, A, B, C,
                                    x_in_idx, y_in_idx, x_edge_idx, y_edge_idx,
                                    x_corner_idx, y_corner_idx,slope_start_idx)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, dt, Lx, Ly, case, True)

        p_saved_l1[:, :, n] = p[x_st_l1_idx: x_en_l1_idx, y_st_l1_idx: y_en_l1_idx]
        p_saved_l2[:, :, n] = p[x_st_l2_idx: x_en_l2_idx, y_st_l1_idx: y_en_l1_idx]

        fsrc1[:, :], fsrc2[:, :] = fsrc.copy(), fsrc1.copy()
        p1[:, :], p2[:, :] = p.copy(), p1.copy()

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'fdtd')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % h_num), t)
    np.save(os.path.join(res_path, 'Ts_%i.npy' % h_num), Ts)
    np.savez_compressed(os.path.join(res_path, 'p_%i_%ideg.npz' % (h_num, int(alpha))),
                        p_l1=p_saved_l1, p_l2=p_saved_l2)
    # np.savez_compressed(os.path.join(res_path, 'p_l1_%i_%ideg.npz' % (h_num, int(alpha))), p_saved_l1)
    # np.savez_compressed(os.path.join(res_path, 'p_l2_%i_%ideg.npz' % (h_num, int(alpha))), p_saved_l2)
