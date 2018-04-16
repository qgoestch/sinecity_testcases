# -*- coding: utf-8 -*-
##
# \file     init_fdtd_paumond_facade.py
# \title    Definition of the numerical parameters for the FDTD method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_fdtd.py files.
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 15 Jan.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fdtd_core')
site.addsitedir(fdtd_core_path)
from upd_fdtd import upd_p_fdtd_srl, upd_vel_pbc_fdtd

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from get_imped_coefts import get_coefts_Miki
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure


def fdtd_srl_init(d_sr, d_sr_ref, h_s, h_r, T, T_delay, f_max_src,
                  rho, c, case, disp_inst_p):
    """
    Setting the 2D geometries and running the FDTD update.
    Main script that contains all the parameters to run the FDTD update in 2D.

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
    :param T_delay: delay of the pulse (s).
    :type T_delay: float
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

    :param      src_typ: source type, string "gauss_1", "sine", "ricker" "dirac"...
    :param      src_frq: source freauency, scalar (Hz).

    :param      dl:      spatial step, scalar (m).
    :param      dt:      time step, scalar (s).
    :param      Lx:      continuous length of the domain (in meter) following the x-direction.
    :param      Ly:      continuous length of the domain (in meter) following the y-direction.
    :param      Nx:      discrete length of the domain (number of node) following the x-direction.
    :param      Ny:      discrete length of the domain (number of node) following the y-direction.
    :param      x:       discrete length sequence of a domain side, scalar (m).
    :param      dx:      spatial step after discretization, scalar (m).
    :param      Nt:      number of iteration time, scalar.
    :param      t:       discretized time sequence, 1d array (s).
    :param      It:      discretized time sequence, 1d array.
    :param      Ts:      time step after dicretization, scalar (s).
    :param      Cn:      Courant number, scalar.

    :param      p:       updated pressure (n+1), numpy array (dimension of the scene).
    :param      p1:      current pressure (n), numpy array (dimension of the scene).
    :param      p2:      past pressure (n-1), numpy array (dimension of the scene).
    :param      fsrc:    soft source (n+1), numpy array (dimension of the scene).
    :param      fsrc1:   soft source (n), numpy array (dimension of the scene).
    :param      p_saved: pressure saved at the receiver location, 1d array (length of the time sequence).

    :param      A:       inertance of the boundary, scalar.
    :param      B:       stiffness of the boundary, scalar.
    :param      C:       resistivity of the boundary, scalar.
    :param      Nb:      boundary of the domain (1 if BC, 0 else) for the compact
                        pressure update, numpy array (dimension of the scene).

    :param      K:       order of the recursive convolution method for the ground impedance, scalar.
    :param      a_k:     residuals of the partial fraction expansion, list of K size.
    :param      gamma_k: poles of the partial fraction expansion, list of K size.

    :param      x_src:   discrete x coordinate of the source, scalar (number of node).
    :param      y_src:   discrete y coordinate of the source, scalar (number of node).
    :param      x_rcv:   discrete x coordinate of the receiver, scalar (number of node).
    :param      y_rcv:   discrete y coordinate of the receiver, scalar (number of node).

    :param      n:       discrete iteration inside the for loop, scalar.

    :return: the acoustic pressure at the pre-defined receivers' locations as a function of time.
    :rtype: (2+1)D array of floats
    """

    # =========================================================================
    #   Size of the numerical domain (mind the reflections !!! NO PML !!!)
    # =========================================================================
    Lx = 30.
    Ly = 5.5 * max(h_s, h_r[-1])

    # =========================================================================
    #   Source
    # =========================================================================
    src_typ = "gauss_1"
    src_frq = f_max_src

    # =========================================================================
    #   Parameters
    # =========================================================================
    dl = 0.075 * c * np.sqrt(2) / src_frq  # 2% accuracy kowalczyk_ieee2011 ~ lambda/9.43
    c = 340.0
    dt = dl / (np.sqrt(2.) * c)

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
    Nt = Nt + nNt   # add the additional time steps to the main time sequence
    print 'Ratio Cn/Cn_th=%g < 1.0' % (Cn * np.sqrt(2))

    src_dly = int(round(T_delay / Ts))

    print '                  FDTD running, details below:                 '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2e s' % Ts
    print 'NUMBER OF It: Nt=%i' % Nt
    print 'DURATION:     T=%.3e s,' % T
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    print 'Sound speed:  c =%.2f m.s-1' % c
    print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m.' % dx
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' \
          % (Nx, Ny, Lx, Ly)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % src_typ
    print 'SOURCE FREQ:  f=%g Hz.' % src_frq
    print 'SOURCE DELAY: %0.2e s, %i n' %(src_dly*Ts, src_dly)
    print '---------------------------------------------------------------'
    # =========================================================================
    #   Variables
    # =========================================================================
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p1 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p2 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc1 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc2 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p_saved_ref = np.zeros((len(d_sr), Nt), dtype=np.complex128)
    p_saved_fac = np.zeros((len(d_sr), len(h_r), Nt), dtype=np.complex128)

    # =========================================================================
    #   Boundaries of the domain, where the BC is calculated
    # =========================================================================
    A = 0.
    B = 1
    C = 0.
    Nb = np.zeros((Nx + 1, Ny + 1))
    i = 1
    Nb[i, 1:-1] = 1.
    i = p.shape[0] - 2
    Nb[i, 1:-1] = 1.
    j = Ny - 1;
    Nb[1:-1, j] = 1.
    j = p.shape[1] - 2
    Nb[1:-1, j] = 1.

    # =========================================================================
    #   Location of the source and receiver(s)
    # =========================================================================
    every_what_pts = 1
    x_src = np.arange(int(round(5. / dx)), int(round(12. / dx)) + 1, every_what_pts)
    x_src1 = int(round(6. / dx))
    x_src2 = int(round(10. / dx))
    y_src = int(round(h_s / dx)) + 2

    x_rcv_0 = [int(round(ii / dx)) for ii in d_sr_ref]
    x_rcv_facade = [int(round(ii / dx)) for ii in d_sr]
    y_rcv = [int(round(ii / dx)) + 2 for ii in h_r]

    # =========================================================================
    #   Calculation of the pressure
    # =========================================================================
    depth = 1
    rand_delay = np.random.randint(len(x_src), size=len(x_src))
    for n in It:
        for ii, x_src_val in enumerate(x_src):
            # fsrc[x_src_val, y_src] = 1. * src.src_select(src_typ, t, n, src_frq,
            #                                             src_dly)
            # fsrc[x_src_val, y_src] = 1. * src.src_select(src_typ, t, n, src_frq,
            #                                       src_dly + int(round(ii * dx * every_what_pts / 340. /Ts)))
            fsrc[x_src_val, y_src] = 1. * src.src_select(src_typ, t, n, src_frq,
                                                  src_dly + int(round(rand_delay[ii] *
                                                              dx * every_what_pts /
                                                              340. /Ts)))
        # fsrc[x_src1, y_src] = 1. * src.src_select(src_typ, t, n, src_frq,
        #                                           src_dly)  # soft source impl.
        # fsrc[x_src2, y_src] = 1. * src.src_select(src_typ, t, n, src_frq,
        #                                           src_dly)
        p = upd_p_fdtd_srl(p, p1, p2, fsrc, fsrc1, fsrc2,
                           Nb, c, rho, Ts, dx, Cn, A, B, C, depth)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, False)

        for d in range(len(d_sr_ref)):
            p_saved_ref[d, n] = p[x_rcv_0[d], y_rcv[0]]

        for d in range(len(d_sr)):
            for hb in range(len(h_r)):
                p_saved_fac[d, hb, n] = p[x_rcv_facade[d], y_rcv[hb]]


        fsrc1[:, :], fsrc2[:, :] = fsrc.copy(), fsrc1.copy()
        p1[:, :], p2[:, :] = p.copy(), p1.copy()

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results',
                            'case%i' % case, 'fdtd')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    np.save(os.path.join(res_path, 't.npy'), t)
    np.save(os.path.join(res_path, 'Ts.npy'), Ts)
    np.save(os.path.join(res_path, 'p_ref.npy'), p_saved_ref)
    np.save(os.path.join(res_path, 'p_fac.npy'), p_saved_fac)
