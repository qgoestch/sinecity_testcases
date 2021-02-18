# -*- coding: utf-8 -*-
##
# \file     init_fdtd_geospr.py
# \title    Definition of the numerical parameters for the FDTD method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_fdtd.py files.
#           It is applied for the grid convergence studies geometrical divergence.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 30 Aug.
##
import numpy as np
import os
import site
import progressbar

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fdtd_core')
site.addsitedir(fdtd_core_path)
from upd_fdtd import upd_p_fdtd_srl, upd_vel_pbc_fdtd

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure


def fdtd_srl_init_conv(dl, h_num, dt, d_sr, T, f_max, rho, c, case, disp_inst_p):
    """
    Setting the 2D geometries and running the FDTD update for case 1: geometrical spreading.
    Main script that contains all the parameters to run the FDTD update in 2D.
    
    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param dt: time step (s).
    :type dt: float
    :param d_sr: horizontal distances between the source and the receivers, list of floats (m).
    :type d_sr: list of floats
    :param T: simulation duration (s).
    :type T: float
    :param f_max: approximated maximale frequency of the source signal (Hz).
    :type f_max: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param disp_inst_p: display the instantaneous pressure.
    :type disp_inst_p: bool
    
    :param Lx: continuous length of the domain (in meter) following the x-direction.
    :param Ly: continuous length of the domain (in meter) following the y-direction.
    :param Nx: discrete length of the domain (number of node) following the x-direction.
    :param Ny: discrete length of the domain (number of node) following the y-direction.
    :param x: discrete length sequence of a domain side, scalar (m).
    :param dx: spatial step after discretization, scalar (m).
    :param Nt: number of iteration time, scalar.
    :param t: discretized time sequence, 1d array (s).
    :param It: discretized time sequence, 1d array.
    :param Ts: time step after dicretization, scalar (s).
    :param Cn: Courant number, scalar.

    :param      p:       updated pressure (n+1), numpy array (dimension of the scene).
    :param      p1:      current pressure (n), numpy array (dimension of the scene).
    :param      p2:      past pressure (n-1), numpy array (dimension of the scene).
    :param      fsrc:    soft source (n+1), numpy array (dimension of the scene).
    :param      fsrc1:   soft source (n), numpy array (dimension of the scene).
    :param      p_axial: pressure saved on the axial axis, 2d array (number of recievers * n).
    :param      p_diag:  pressure saved on the diagonal axis, 2d array (number of recievers * n).

    :param      A:       inertance of the boundary, scalar.
    :param      B:       stiffness of the boundary, scalar.
    :param      C:       resistivity of the boundary, scalar.
    :param      Nb:      boundary of the domain (1 if BC, 0 else) for the compact
                        pressure update, numpy array (dimension of the scene).

    :param      x_src:   discrete x coordinate of the source, scalar (number of node).
    :param      y_src:   discrete y coordinate of the source, scalar (number of node).
    :param      x_rcv:   discrete x coordinate of the receiver, scalar (number of node).
    :param      y_rcv:   discrete y coordinate of the receiver, scalar (number of node).

    :param      n:       discrete iteration inside the for loop, scalar.

    :return: the acoustic pressure at the pre-defined receivers' locations as a function of time.
    :rtype: (1+1)D array of floats 64
    """

    # ==============================================================================
    #   Source
    # ==============================================================================
    src_typ = "gauss_1"
    src_frq = f_max

    # ==============================================================================
    #   Parameters
    # ==============================================================================
    dl_max = 0.075 * c * np.sqrt(2) / src_frq  # 2% accuracy kowalczyk_ieee2011 ~ lambda/9.43
    dt = dl / (np.sqrt(2.) * c)
    # print 'dt_lim = %.3e' %dt_lim
    Lx = 24.0
    Ly = Lx # compulsory for the right diagonal coordinates
    Nx = np.int(np.round(Lx / dl))
    Ny = np.int(np.round(Ly / dl))
    x = np.linspace(0, Lx, Nx + 1)
    dx = np.float64(x[1] - x[0])
    dx = round(dx, 5)
    Nt = int(round(T / float(dt)))
    t = np.linspace(0, Nt * dt, Nt + 1)
    It = range(0, t.shape[0])
    Ts = np.float64(t[1] - t[0])
    Cn = np.float64(c * Ts / dx)
    nNt = 1
    # print '--------------------- Courant Number --------------------------'
    # while Cn >= (1 / np.sqrt(2)):
    #     print 'COURANT NUMBER CORRECTION!!!'
    #     nNt = 1 + nNt
    #     t = np.linspace(0, Nt * dt, Nt + (nNt))  # time discretization
    #     It = range(0, t.shape[0] - 1)  # time iterations range
    #     Ts = np.float64(t[1] - t[0])  # sampling period for staggered grid
    #     Cn = np.float64(c * Ts / dx)
    # print 'Additional iterations for stab: %i' % (nNt)
    # Nt = Nt + nNt  # add the additional time steps to the main time sequence
    # print 'Ratio Cn/Cn_th=%g < 1.0' % (Cn * np.sqrt(2))

    src_dly = int(round(1/50./Ts))

    print '                 FDTD geometrical spreading                    '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2g s' % (Ts)
    print 'NUMBER OF It: Nt=%i' % (Nt)
    print 'DURATION:     T=%.3f s,' % (T)
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    print '2PERCENT ACCURACY:  %.3f Hz,' % (0.075 / Ts)
    print 'COURANT NUMBER: Cn = %.3f' % Cn
    print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m, dl_max=%g m.' % (dx, dl_max)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' \
          % (Nx, Ny, Lx, Ly)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % src_typ
    print 'SOURCE FREQ:  f=%g Hz.' % src_frq
    print 'SOURCE DELAY:  i=%g Iter.' % src_dly
    print '---------------------------------------------------------------'
    # ==============================================================================
    #   Variables
    # ==============================================================================
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p1 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p2 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc1 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    fsrc2 = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p_axial = np.zeros((len(d_sr), Nt + 1), dtype=np.complex128)
    p_diag  = np.zeros((len(d_sr), Nt + 1), dtype=np.complex128)

    # ==============================================================================
    #   Boundaries of the domain, where the BC is calculated
    # ==============================================================================
    A = 0.;
    B = 1;
    C = 0.
    Nb = np.zeros((Nx + 1, Ny + 1))
    i = 1;
    Nb[i, 1:-1] = 1.
    i = p.shape[0] - 2;
    Nb[i, 1:-1] = 1.
    j = Ny - 1;
    Nb[1:-1, j] = 1.
    j = p.shape[1] - 2;
    Nb[1:-1, j] = 1.

    # ==============================================================================
    #   Location of the source and receiver(s)
    # ==============================================================================
    x_src = int(round(Lx / 2. / dx))
    y_src = int(round(Ly / 2. / dx))

    x_axial = [int(round((Lx / 2. + ii) / dx )) for ii in d_sr]
    y_axial = y_src

    x_diag = [int(round( (Lx / 2. + ii*np.cos(np.pi/4)) / dx )) for ii in d_sr]
    y_diag = x_diag

    # # ==============================================================================
    # #   Remove the old figures of wavefronts
    # # ==============================================================================
    # if disp_inst_p:
    #     import os
    #     import glob
    #     res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case%i' % case, 'figures')
    #     for name in glob.glob(os.path.join(res_path, 'wave_front_*')):
    #         os.remove(name)  # clean up previous plot files

    # ==============================================================================
    #   Calculation of the pressure
    # ==============================================================================
    depth = 1
    bar = progressbar.ProgressBar(maxval=It[-1],
                                  widgets={'(', progressbar.Timer(), ') ', progressbar.Bar('x'), ' (', progressbar.ETA(), ') ', })
    bar.start()
    for n in It:
        # print("Iteration rate %i/%i" % (n, max(It)))
        # p1[x_src, y_src] = 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # hard source impl.
        fsrc[x_src, y_src] = 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # soft source impl.
        p = upd_p_fdtd_srl(p, p1, p2, fsrc, fsrc1, fsrc2,
                           Nb, c, rho, Ts, dx, Cn, A, B, C, depth)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, True)

        for d in range(len(d_sr)):
            p_axial[d, n] = p[x_axial[d], y_axial]
            p_diag[d, n] = p[x_diag[d], y_diag[d]]

        fsrc1[:, :], fsrc2[:, :] = fsrc.copy(), fsrc1.copy()
        p1[:, :], p2[:, :] = p.copy(), p1.copy()
        bar.update(n)
    bar.finish()
    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'fdtd')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % (h_num)), t)
    np.save(os.path.join(res_path, 'Ts_%i.npy' % (h_num)), Ts)
    np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('axial', h_num)), p_axial)
    np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('diag', h_num)), p_diag)
