# -*- coding: utf-8 -*-
##
# \file     init_tlm_modes.py
# \title    Definition of the numerical parameters for the TLM method for the
#           grid convergence sudy. The satial is an input parameter.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_tlm.py files.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 26 Jul.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

tlm_core_path = os.path.join(base_path, 'tlm_core')
site.addsitedir(tlm_core_path)
from upd_tlm import upd_p_tlm_srl_rigid

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_modes

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import inst_pres_exact_vs_num


def tlm_srl_init_conv_box(dl, h_num, h_set, t_save, T, Lx, Ly, nx, ny, fr, rho, c, case, disp_inst_p):
    """

         Setting the 2D geometries and running the TLM update for case 2: acoustic modes.
      Main script that contains all the parameters to run the TLM update in 2D.

    :param dl: spatial step (m).
    :type dl: float
    :param h_num: spatial step index.
    :type h_num: int
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param t_save: time instant at which the pressure maps are saved for comparison (s).
    :type t_save: float
    :param T: simulation duration (s).
    :type T: float
    :param Lx: length of the box followinf the x axis (m)
    :type Lx: float
    :param Ly: length of the box followinf the y axis (m)
    :type Ly: float
    :param nx: mode number following the x direction
    :type nx: int
    :param ny: mode number following the y direction
    :type ny: int
    :param fr:
    :type fr:
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param disp_inst_p: display the instantaneous pressure.
    :type disp_inst_p: bool

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

    :param      I_t     incident pulse from the top.
    :param      I_b     incident pulse from the bottom.
    :param      I_l     incident pulse from the left.
    :param      I_r     incident pulse from the right.

    :param      geo     boundary of the domain (1 if BC, 0 else), similar to Nb
                        in FDTD scripts, numpy array (dimension of the scene).

    :return: the acoustic pressure at the pre-defined receivers' locations as a function of time.
    :rtype: (1+1)D array of floats 64
    """
    # ==============================================================================
    #   Parameters
    # ==============================================================================
    # dt = dl / (np.sqrt(2.) * c)
    Cn_lim = np.sqrt(2)**-1
    dt_coarse = 2 * 10 ** -4
    dl_coarse = h_set[-1]
    c = np.float64(Cn_lim * dl_coarse / dt_coarse)
    print 'sound speed: c=%f m.s-1' % c
    dt = np.float64(Cn_lim * dl / c)

    Lx = Lx + dl
    Ly = Ly + dl
    Nx = np.int(np.round(Lx / dl))
    Ny = np.int(np.round(Ly / dl))
    x = np.linspace(0, Lx, Nx + 1)
    y = np.linspace(0, Ly, Ny + 1)
    dx = np.float64(x[1] - x[0])
    # dx = round(dx, 5)
    Nt = int(round(T / float(dt)))
    t = np.linspace(0, Nt * dt, Nt + 1)
    It = range(0, t.shape[0])
    Ts = np.float64(t[1] - t[0])


    print '                     TLM in 2D box                        '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.3g s' % (Ts)
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m.' % (dx)
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' % (Nx, Ny, Lx, Ly)

    # ==============================================================================
    #   Incident (I) and Reflected (R) pulses from all directions (top, bottom...)
    # ==============================================================================
    I_t = np.zeros((Nx + 1, Ny + 1))
    I_b = np.zeros((Nx + 1, Ny + 1))
    I_l = np.zeros((Nx + 1, Ny + 1))
    I_r = np.zeros((Nx + 1, Ny + 1))

    p = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_saved = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_exact = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_exact_saved = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)

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
    #   Definition of the initial conditions
    # ==============================================================================
    A_tlm = 0.5
    p_exact[1:-1, 1:-1] = analytic_solution_modes(dx,p,nx,ny,x,y,c,t,0)
    I_t[1:-1, 1:-1] = A_tlm * p_exact[1:-1, 1:-1]
    I_b[1:-1, 1:-1] = A_tlm * p_exact[1:-1, 1:-1]
    I_l[1:-1, 1:-1] = A_tlm * p_exact[1:-1, 1:-1]
    I_r[1:-1, 1:-1] = A_tlm * p_exact[1:-1, 1:-1]
    # it_save = int(round(4*dt_coarse/ Ts))
    it_save = int(round(t_save / Ts))
    print 't_save=%0.3e s; it_save = %i it' %(4*dt_coarse,it_save)
    for n in It[:-1]:
        p_exact[1:-1, 1:-1] = analytic_solution_modes(dx,p,nx,ny,x,y,c,t,n)
        # ==============================================================================
        #   Calculation of the Incident and Scattered pulses (via upd_tlm.py)
        # ==============================================================================
        p = upd_p_tlm_srl_rigid(I_t,I_b,I_l,I_r)

        if disp_inst_p:
            inst_pres_exact_vs_num(n, p_exact, p, Lx, Ly, nx, ny, dx, dt)

        if n == it_save + 0:
            p_exact_saved[1:-1,1:-1] = p_exact[1:-1,1:-1]

        if n == it_save:
            p_saved[1:-1,1:-1] = p[1:-1,1:-1]

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i' % case, 'tlm')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % (h_num)), t)
    # np.save(os.path.join(res_path, 'Ts_%i.npy' % (h_num)), Ts)
    np.save(os.path.join(res_path, 'p_tlm_%i.npy' % (h_num)), p_saved)
    np.save(os.path.join(res_path, 'p_an_%i.npy' % (h_num)), p_exact_saved)