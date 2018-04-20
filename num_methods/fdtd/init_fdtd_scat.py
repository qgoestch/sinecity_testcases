# -*- coding: utf-8 -*-
##
# \file     init_fdtd_scat.py
# \title    Definition of the numerical parameters for the FDTD method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_fdtd.py files.
#           It is applied for the grid convergence studies of 2D scattering by a circular obstacle.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 07 Sep.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fdtd_core')
site.addsitedir(fdtd_core_path)
from upd_fdtd import upd_p_fdtd_srl, upd_p_fdtd_srl_2D_scat

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
import source_signals as src
from D_circ_final import D_circ_surf

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure
from plot_boundary_cells import plot_cells_circle

def fdtd_srl_init_scat(dl, h_num, radius, T_delay, T, f_max, rho, c, case, disp_inst_p,free_field):
    """
         Setting the 2D geometries and running the FDTD update for case 4: scattering.
      Main script that contains all the parameters to run the FDTD update in 2D.

    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param radius: radius of the scatterer (m)
    :type radius: float
    :param T_delay: delay of the pulse (s)
    :type T_delay: float
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
    :param free_field: the domain is enlarged
    :type free_field: bool

    :param      src_typ source type, string "gauss_1", "sine", "ricker" "dirac"...
    :param      src_frq source freauency, scalar (Hz).

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

    :param      x_in_idx    x coordinates of the zero cells inside the cylinder, list of integers.
    :param      y_in_idx    y coordinates of the zero cells inside the cylinder, list of integers.
    :param      x_edges_idx x coordinates of the edge cells with 1 faces connected to the circle, list of integers.
    :param      y_edges_idx y coordinates of the edge cells with 1 faces connected to the circle, list of integers.
    :param      x_corners_idx   x coordinates of the edge cells with 2 faces connected to the circle, list of integers.
    :param      y_corners_idx   y coordinates of the edge cells with 2 faces connected to the circle, list of integers.

    :param      p       updated pressure (n+1), numpy array (dimension of the scene).
    :param      p1      current pressure (n), numpy array (dimension of the scene).
    :param      p2      past pressure (n-1), numpy array (dimension of the scene).
    :param      fsrc    soft source (n+1), numpy array (dimension of the scene).
    :param      fsrc1   soft source (n), numpy array (dimension of the scene).

    :param      A       inertance of the boundary, scalar.
    :param      B       stiffness of the boundary, scalar.
    :param      C       resistivity of the boundary, scalar.
    :param      Nb      boundary of the domain (1 if BC, 0 else) for the compact
                        pressure update, numpy array (dimension of the scene).

    :param      rad_rcv radius of the receivers circles, float (m).
    :param      phi_rcv angles of each discrete receivers on the circles, float (rad).
    :param      ipm     x discrete corrdinates -grid index- of the receivers, list of integers.
    :param      jpm     y discrete corrdinates -grid index- of the receivers, list of integers.
    :param      rcv_dist_idx    exact distance on the grid, used for the analytic solution, float (m).
    :param      rcv_phi_idx     exact angles on the grid, used for the analytic solution, float (rad).
    :param      p_circ  pressure saved at the receiver location, on circles, 3d array (distance, angles, time).

    :param      n       discrete iteration inside the for loop, scalar.

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
    Lx = 12.
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
    nNt = 0
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

    src_dly = int(T_delay/Ts)

    print '                     FDTD 2D circular obstacle                 '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2g s' % (Ts)
    print 'NUMBER OF It: Nt=%i' % (Nt)
    print 'DURATION:     T=%.3f s,' % (T)
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    print '2PERCENT ACCURACY:  %.3f Hz,' % (0.075 / Ts)
    print 'COURANT NUMBER: Cn = %.3f' % (Cn)
    print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m, dl_max=%g m.' % (dx, dl_max)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' % (Nx, Ny, Lx, Ly)
    #    print 'IMPEDANCE:    Yb=%f .' %(Yb)
    print '---------------------- Source Signal --------------------------'
    print 'SOURCE TYPE:  %s,' % (src_typ)
    print 'SOURCE FREQ:  f=%g Hz.' % (src_frq)
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

    # ==============================================================================
    #   Definition of the circle
    # ==============================================================================
    ray = radius
    if not free_field:
        x_corners_idx,      y_corners_idx,\
        x_edges_idx,        y_edges_idx, \
        x_in_idx,           y_in_idx, \
        x_corners_idx_1st, y_corners_idx_1st, x_corners_idx_2nd, y_corners_idx_2nd, \
        x_corners_idx_3rd, y_corners_idx_3rd, x_corners_idx_4th, y_corners_idx_4th, \
        x_edges_idx_1st, y_edges_idx_1st, x_edges_idx_2nd, y_edges_idx_2nd, \
        x_edges_idx_3rd, y_edges_idx_3rd, x_edges_idx_4th, y_edges_idx_4th \
        = D_circ_surf(dx,ray,Lx,Ly)

        bc_ed = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
        bc_co = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
        bc_in = np.zeros((Nx + 1, Ny + 1), dtype=np.int64)
        for ndx in range(len(x_in_idx)):
            bc_in[x_in_idx[ndx], y_in_idx[ndx]] = 4
        for ndx in range(len(x_edges_idx)):
            bc_ed[x_edges_idx[ndx], y_edges_idx[ndx]] = 1
        for ndx in range(len(x_corners_idx)):
            bc_co[x_corners_idx[ndx], y_corners_idx[ndx]] = 2
        # plot_cells_circle(dx, bc_ed, bc_co, bc_in, ray, Nx, Ny)

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
    #   Location of the receivers
    # ==============================================================================
    ax = np.floor(Nx / 2) + dx
    by = np.floor(Ny / 2)
    rad_rcv = np.linspace(0.5, 3.0, 26)
    phi_rcv = np.linspace(0, 2 * np.pi, 120)
    ipm = np.zeros((len(phi_rcv), len(rad_rcv)), dtype=np.int)
    jpm = np.zeros((len(phi_rcv), len(rad_rcv)), dtype=np.int)
    rcv_dist_idx = np.zeros((len(phi_rcv), len(rad_rcv)))
    rcv_phi_idx = np.zeros((len(phi_rcv), len(rad_rcv)))
    for iangle in range(len(phi_rcv)):
        for idist in range(len(rad_rcv)):
            ipm[iangle, idist] = int(round(ax + np.around(rad_rcv[idist] / dx * np.cos(phi_rcv[iangle]))))
            jpm[iangle, idist] = int(round(by + np.around(rad_rcv[idist] / dx * np.sin(phi_rcv[iangle]))))
            rcv_dist_idx[iangle, idist] = np.sqrt(
                (np.abs(ipm[iangle, idist] - ax) * dx) ** 2 +
                (np.abs(jpm[iangle, idist] - by) * dx) ** 2)
            if iangle <= int(len(phi_rcv) / 2):
                #                print 'yes'
                rcv_phi_idx[iangle, idist] = np.abs(np.arccos((ax - ipm[iangle, idist]) * dx /
                                                               rcv_dist_idx[iangle, idist]) - np.pi)

            elif iangle > int(len(phi_rcv) / 2):
                #                print 'no'
                rcv_phi_idx[iangle, idist] = np.arccos((ax - ipm[iangle, idist]) * dx /
                                                        rcv_dist_idx[iangle, idist]) + np.pi
    p_circ = np.zeros((len(phi_rcv), len(rad_rcv), Nt + 1), dtype=np.complex128)

    # ==============================================================================
    #   Calculation of the pressure
    # ==============================================================================
    for n in It:
        p1[1,1:-1]= 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # hard source impl.
        # fsrc[1,1:-1] = 1. * src.src_select(src_typ, t, n, src_frq, src_dly)  # soft source impl.

        if not free_field:
            p = upd_p_fdtd_srl_2D_scat( p, p1, p2, fsrc, fsrc2,
                                        Nb, c, rho, Ts, dx, Cn, A, B, C,
                                        x_in_idx, y_in_idx, x_edges_idx, y_edges_idx,
                                        x_corners_idx, y_corners_idx)
        else:
            p = upd_p_fdtd_srl( p,p1,p2,fsrc,fsrc1,fsrc2,Nb,c,rho,Ts,dx,Cn,A,B,C,1)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, True)

        for m_ang in range(len(phi_rcv)):
            for n_dist in range(len(rad_rcv)):
                p_circ[m_ang,n_dist,n]=p[ipm[m_ang,n_dist],jpm[m_ang,n_dist]]

        fsrc1[:, :], fsrc2[:, :] = fsrc.copy(), fsrc1.copy()
        p1[:, :], p2[:, :] = p.copy(), p1.copy()

    # ==============================================================================
    #   Save the variables used for the calculation of the error
    # ==============================================================================
    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'fdtd','p_h%i' %(h_num))
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    if free_field:
        field = 'f'
        np.save(os.path.join(res_path, 'Ts.npy'  ), Ts)
        np.save(os.path.join(res_path, 'h.npy'   ), dx)
        np.save(os.path.join(res_path, 't.npy'), t)
        np.save(os.path.join(res_path, 'p_%s.npy'%(field)), p_circ)
    else:
        field = 't'
        np.save(os.path.join(res_path, 'rcvdist.npy' ), rcv_dist_idx)
        np.save(os.path.join(res_path, 'rcvphi.npy'  ), rcv_phi_idx)
        np.save(os.path.join(res_path, 'p_%s.npy' %(field)), p_circ)