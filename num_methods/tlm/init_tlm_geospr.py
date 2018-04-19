# -*- coding: utf-8 -*-
##
# \file     init_tlm_geospr.py
# \title    Definition of the numerical parameters for the TLM method for the
#           grid convergence sudy. The satial is an input parameter.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_tlm.py files.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 01 Sep.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

tlm_core_path = os.path.join(base_path, 'tlm_core')
site.addsitedir(tlm_core_path)
from upd_tlm import upd_p_tlm_srl_rigid

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import instatenous_pressure

def tlm_srl_init_conv(dl, h_num, dt, d_sr, T, f_max, rho, c, case, disp_inst_p):
    """
         Setting the 2D geometries and running the TLM update for case 1: geometrical spreading.
      Main script that contains all the parameters to run the TLM update in 2D.

    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param dt: time step (s).
    :type dt: float
    :param d_sr: horizontal distances between the source and the receivers (m).
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

    :param      x_src   discrete x coordinate of the source, scalar (number of node).
    :param      y_src   discrete y coordinate of the source, scalar (number of node).
    :param      x_rcv   discrete x coordinate of the receiver, scalar (number of node).
    :param      y_rcv   discrete y coordinate of the receiver, scalar (number of node).

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

    :param      p_axial pressure saved on the axial axis, 2d array (number of recievers * n).
    :param      p_diag  pressure saved on the diagonal axis, 2d array (number of recievers * n).

    :return: the acoustic pressure at the pre-defined receivers' locations as a function of time.
    :rtype: (2+1)D array of floats
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

    src_dly = int(round(1 / 50. / Ts))

    # ==============================================================================
    #   Location of the source and receiver(s)
    # ==============================================================================
    x_src = int(round(Lx / 2. / dx))
    y_src = int(round(Ly / 2. / dx))

    x_axial = [int(round((Lx / 2. + ii) / dx )) for ii in d_sr]
    y_axial = y_src

    x_diag = [int(round( (Lx / 2. + ii*np.cos(np.pi/4)) / dx )) for ii in d_sr]
    y_diag = x_diag

    print '                  TLM geometrical spreading                    '
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
    S_t = np.zeros((Nx + 1, Ny + 1))
    S_b = np.zeros((Nx + 1, Ny + 1))
    S_l = np.zeros((Nx + 1, Ny + 1))
    S_r = np.zeros((Nx + 1, Ny + 1))
    p = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    p_axial = np.zeros((len(d_sr), Nt + 1), dtype=np.complex128)
    p_diag  = np.zeros((len(d_sr), Nt + 1), dtype=np.complex128)

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
    # Normalization of the amplitude according to the FDTD signal at 0.5 m, see plot_time_signals.py (last figure)
    # A_tlm = 0.5 * 1. # 500 Hz
    A_tlm = 0.5 * 5.  # 2000 Hz
    for n in It[:-1]:
        # source signal impl.
        I_t[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_b[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_l[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)
        I_r[x_src, y_src] = A_tlm * src.src_select(src_typ, t, n, src_frq, src_dly)

        # ==============================================================================
        #   Calculation of the Incident and Scattered pulses (via upd_tlm.py)
        # ==============================================================================
        p = upd_p_tlm_srl_rigid(I_t, I_b, I_l, I_r)

        if disp_inst_p:
            instatenous_pressure(n, Nt, p, dx, Ts, Lx, Ly, case, False)

        for d in range(len(d_sr)):
            p_axial[d, n] = p[x_axial[d], y_axial]
            p_diag[d, n] = p[x_diag[d], y_diag[d]]

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'tlm')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % (h_num)), t)
    np.save(os.path.join(res_path, 'Ts_%i.npy' % (h_num)), Ts)
    np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('axial', h_num)), p_axial)
    np.save(os.path.join(res_path, 'p_%s_%i.npy' % ('diag', h_num)), p_diag)