# -*- coding: utf-8 -*-
##
# \file     init_fdtd_modes.py
# \title    Definition of the numerical parameters for the FDTD method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_fdtd.py files.
#           It is applied for the grid convergence studies for modes of a 2D square box.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 25 Jul.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_core_path = os.path.join(base_path, 'fdtd_core')
site.addsitedir(fdtd_core_path)
from upd_fdtd import upd_p_fdtd_srl, upd_vel_pbc_fdtd

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_modes

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from display_wavefronts import inst_pres_exact_vs_num


def fdtd_srl_init_conv_box(dl, h_num, h_set, t_save, T, Lx, Ly, nx, ny, fr, rho, c, case, disp_inst_p):
    """

         Setting the 2D geometries and running the FDTD update for case 2: acoustic modes.
      Main script that contains all the parameters to run the FDTD update in 2D.

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

    :param      n       discrete iteration inside the for loop, scalar.

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
    dx = round(dx, 5)
    Nt  = int(round(T/dt))                                                      # number of iterations
    t   = np.linspace(0, Nt*dt, Nt+1)                                # time discretization
    It  = range(len(t))                                                         # time grid
    Ts  = np.float64(t[1] - t[0])                                                # time step (from the time )
    Cn  = np.float64((c*Ts/dx))
    # nNt = 0
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
    print 'Ratio Cn/Cn_th=%g < 1.0' % (Cn * np.sqrt(2))

    print '                     FDTD in 2D box                        '
    print '-------------------------- Time -------------------------------'
    print 'TIME-STEP:    Ts=%0.3g s' % Ts
    print 'NUMBER OF It: Nt=%i' % Nt
    print 'DURATION:     T=%.3f s,' % T
    print '-------------------------- Space ------------------------------'
    print 'SPATIAL-STEP: dx=%g m.' % dx
    print 'DIMENSIONS:   Nx=%i cells; Ny=%i cells; Lx=%g m; Ly=%g m.' \
          % (Nx, Ny, Lx, Ly)

    # =========================================================================
    #   Variables
    # =========================================================================
    p       = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p1      = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p2      = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    fsrc    = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    fsrc1   = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    fsrc2   = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_saved = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_exact = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    p_exact_saved = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)


    p_bc    = np.zeros((Nx + 1, Ny + 1))
    p1_bc   = np.zeros((Nx + 1, Ny + 1))

    v_x     = np.zeros((Nx + 1, Ny + 1))
    v1_x    = np.zeros((Nx + 1, Ny + 1))
    v_y     = np.zeros((Nx + 1, Ny + 1))
    v1_y    = np.zeros((Nx + 1, Ny + 1))

    # =========================================================================
    #   Boundaries of the domain, where the BC is calculated
    # =========================================================================
    A = 0.;
    B = 1;
    C = 0.
    Nb = np.zeros((Nx + 1, Ny + 1))
    BC_node = 1
    """ Edges index """
    i=BC_node;                  Nb[i,1:-1]=1.
    i=p.shape[0]-BC_node-1;     Nb[i,1:-1]=1.
    j=BC_node;                  Nb[1:-1,j]=1.
    j=p.shape[1]-BC_node-1;     Nb[1:-1,j]=1.

    """ Corner index """
    for i in [BC_node,p.shape[0]-BC_node-1]:
        for j in [BC_node,p.shape[1]-BC_node-1]:
            Nb[i,j]=2.

    # =========================================================================
    #   Definition of the initial conditions
    # =========================================================================
    p2[1:-1,1:-1] = analytic_solution_modes(dx,p,nx,ny,x,y,c,t,0)
    p1[1:-1,1:-1] = analytic_solution_modes(dx,p,nx,ny,x,y,c,t,1)
    p = upd_p_fdtd_srl(p,p1,p2,fsrc,fsrc1,fsrc2,Nb,c,rho,Ts,dx,Cn,A,B,C,1)
    p1[:,:], p2[:,:] = p.copy(), p1.copy()

    # =========================================================================
    #   Calculation of the pressure
    # =========================================================================
    depth = 1
    # it_save = int(round(4*dt_coarse/ Ts))
    it_save = int(round(t_save / Ts))
    print 't_save=%0.3e s; it_save = %i it' %(4*dt_coarse,it_save)
    for n in It[2:-1]:
        p_exact[1:-1,1:-1] = analytic_solution_modes(dx,p,nx,ny,x,y,c,t,n)
        p = upd_p_fdtd_srl(p, p1, p2, fsrc, fsrc1, fsrc2,
                           Nb, c, rho, Ts, dx, Cn, A, B, C, depth)

        if disp_inst_p:
            inst_pres_exact_vs_num(n, p_exact, p, Lx, Ly, nx, ny, dx, dt)

        if n == it_save + 0:
            print '----> Saving exact <----- t[n]=%f, It=%i' %(t[n],n)
            print np.max(p_exact[:,:])#,np.min(p_exact[:,:,:]),np.average(p_exact[:,:,:])
            p_exact_saved[1:-1,1:-1] = p_exact[1:-1,1:-1]

        if n == it_save:
            print '----> Saving FDTD <----- t[n]=%f, It=%i' %(t[n],n)
            print np.max(p[:,:])#,np.min(p_exact[:,:,:]),np.average(p_exact[:,:,:])
            p_saved[1:-1,1:-1] = p[1:-1,1:-1]

        fsrc1[:, :], fsrc2[:, :] = fsrc.copy(), fsrc1.copy()
        p1[:, :], p2[:, :] = p.copy(), p1.copy()
        p1_bc[:, :], v1_x[:, :], v1_y[:, :] = p_bc.copy(), v_x.copy(), v_y.copy()

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'results', 'case%i'
                            % case, 'fdtd')
    if not os.path.exists(res_path): os.makedirs(res_path)
    np.save(os.path.join(res_path, 't_%i.npy' % (h_num)), t)
    # np.save(os.path.join(res_path, 'Ts_%i.npy' % (h_num)), Ts)
    np.save(os.path.join(res_path, 'p_fdtd_%i.npy' % (h_num)), p_saved)
    np.save(os.path.join(res_path, 'p_exact_%i.npy' % (h_num)), p_exact)
    np.save(os.path.join(res_path, 'p_an_%i.npy' % (h_num)), p_exact_saved)