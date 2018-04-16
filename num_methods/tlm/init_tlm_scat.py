# -*- coding: utf-8 -*-
##
# \file     init_tlm_scat.py
# \title    Definition of the numerical parameters for the TLM method for the
#           grid convergence sudy. The satial is an input parameter.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_tlm.py files.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 08 Sep.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

tlm_core_path = os.path.join(base_path, 'tlm_core')
site.addsitedir(tlm_core_path)
from upd_tlm import upd_p_tlm_srl_rigid, upd_p_tlm_srl_2D_scat

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from get_imped_coefts import get_coefts_Miki
import source_signals as src
from D_circ_final import D_circ_surf

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure

def tlm_srl_init_scat(dl, h_num, radius, T_delay, T, f_max, rho, c, case, disp_inst_p,free_field):
    """
         Setting the 2D geometries and running the TLM update for case 4: scattering.
      Main script that contains all the parameters to run the TLM update in 2D.

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

    :param      rad_rcv radius of the receivers circles, float (m).
    :param      phi_rcv angles of each discrete receivers on the circles, float (rad).
    :param      ipm     x discrete corrdinates -grid index- of the receivers, list of integers.
    :param      jpm     y discrete corrdinates -grid index- of the receivers, list of integers.
    :param      rcv_dist_idx    exact distance on the grid, used for the analytic solution, float (m).
    :param      rcv_phi_idx     exact angles on the grid, used for the analytic solution, float (rad).
    :param      p_circ  pressure saved at the receiver location, on circles, 3d array (distance, angles, time).

    :param      I_t     incident pulse from the top.
    :param      I_b     incident pulse from the bottom.
    :param      I_l     incident pulse from the left.
    :param      I_r     incident pulse from the right.
    :param      S_t     scattered pulse from the top.
    :param      S_b     scattered pulse from the bottom.
    :param      S_l     scattered pulse from the left.
    :param      S_r     scattered pulse from the right.

    :param      geo     boundary of the domain (1 if BC, 0 else), similar to Nb
                        in FDTD scripts, numpy array (dimension of the scene).

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

    src_dly = int(T_delay/Ts)

    # ==============================================================================
    #   Definition of the circle
    # ==============================================================================
    ray = radius
    if not free_field:
        x_corners_idx, y_corners_idx, \
        x_edges_idx, y_edges_idx, \
        x_in_idx, y_in_idx, \
        x_corners_idx_1st, y_corners_idx_1st, x_corners_idx_2nd, y_corners_idx_2nd, \
        x_corners_idx_3rd, y_corners_idx_3rd, x_corners_idx_4th, y_corners_idx_4th, \
        x_edges_idx_1st, y_edges_idx_1st, x_edges_idx_2nd, y_edges_idx_2nd, \
        x_edges_idx_3rd, y_edges_idx_3rd, x_edges_idx_4th, y_edges_idx_4th \
        = D_circ_surf(dx,ray,Lx,Ly)

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

    print '                     TLM 2D circular obstacle                  '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.2g s' % (Ts)
    print 'NUMBER OF It: Nt=%i' % (Nt)
    print 'DURATION:     T=%.3f s,' % (T)
    print 'SAMP. FREQ.:  Fs=%.3f Hz,' % (1 / Ts)
    print 'BANDWIDTH:    FMAX=%.3f Hz,' % (0.196 / Ts)
    print '2PERCENT ACCURACY:  %.3f Hz,' % (0.075 / Ts)
    print 'Nodes/lambda = %.3f' % (c / src_frq / dx)
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m, dl_max=%g m.' % (dx, dl_max)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' % (Nx, Ny, Lx, Ly)
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
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)

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
    #   Hard source assignment
    # ==============================================================================
    A_fdtd_inv = 1.e0  # see gauss_1 in source_signals.py
    A_tlm = 0.5 * A_fdtd_inv
    for n in It[:-1]:
        # source signal impl.
        I_t[1,1:-1] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_b[1,1:-1] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_l[1,1:-1] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_r[1,1:-1] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)

        # ==============================================================================
        #   Calculation of the Incident and Scattered pulses (via upd_tlm.py)
        # ==============================================================================
        if not free_field:
            p = upd_p_tlm_srl_2D_scat(I_t, I_b, I_l, I_r, x_in_idx, y_in_idx,
                                      x_edges_idx_1st, y_edges_idx_1st, x_edges_idx_2nd, y_edges_idx_2nd,
                                      x_edges_idx_3rd, y_edges_idx_3rd, x_edges_idx_4th, y_edges_idx_4th,
                                      x_corners_idx_1st, y_corners_idx_1st, x_corners_idx_2nd, y_corners_idx_2nd,
                                      x_corners_idx_3rd, y_corners_idx_3rd, x_corners_idx_4th, y_corners_idx_4th)
        else:
            p = upd_p_tlm_srl_rigid(I_t, I_b, I_l, I_r)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, False)

        for m_ang in range(len(phi_rcv)):
            for n_dist in range(len(rad_rcv)):
                p_circ[m_ang,n_dist,n]=p[ipm[m_ang,n_dist],jpm[m_ang,n_dist]]

    # ==============================================================================
    #   Save the variables used for the calculation of the error
    # ==============================================================================
    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'tlm','p_h%i' %(h_num))
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    if free_field:
        field = 'f'
        np.save(os.path.join(res_path, 'Ts.npy'  ), Ts)
        np.save(os.path.join(res_path, 'h.npy'   ), dx)
        np.save(os.path.join(res_path, 't.npy' ), t)
        np.save(os.path.join(res_path, 'p_%s.npy'%(field)), p_circ)
    else:
        field = 't'
        np.save(os.path.join(res_path, 'rcvdist.npy' ), rcv_dist_idx)
        np.save(os.path.join(res_path, 'rcvphi.npy'  ), rcv_phi_idx)
        np.save(os.path.join(res_path, 'p_%s.npy' %(field)), p_circ)