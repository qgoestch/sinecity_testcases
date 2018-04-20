# -*- coding: utf-8 -*-
##
# \file     display_wavefronts.py
# \title    Comparison of the spectrum and instantaneous pressure display.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 06 Oct.
##
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
base_path = reduce (lambda l,r: l + os.path.sep + r, 
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )

def inst_pres_exact_vs_num(n,p_exact,p_num, Lx, Ly, nx, ny, dx, dt):
    """
    Display the instaneous pressure in init_*.py if disp_inst_p == True.

    :param n: time iteration index
    :type n: int
    :param p_exact: analytic solution for the acoustic pressure (Pa)
    :type p_exact: 2D numpy array of floats.
    :param p_num: numerical results for the acoustic pressure (Pa)
    :type p_num: 2D numpy array of floats.
    :param Lx: continuous length of the domain (in meter) following the x-direction
    :type Lx: float
    :param Ly: continuous length of the domain (in meter) following the y-direction
    :type Ly: float
    :param nx: mode number for the x direction
    :type nx: int
    :param ny: mode number for the y direction
    :type ny: int
    :param dx: spatial step (m).
    :type dx: float
    :param dt: time step (s).
    :type dt: float
    :return: display (and save snapchots of) the instantaneous pressure in the 2D domain as a function of time.
    """
    if n > 1 and n % 5 == 0:  # instant pressure plot
        fig = plt.figure(1)
        plt.clf()
        ax1 = plt.subplot(211)
        h_x = 40
        h_y = h_x
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(h_x))
        ax1.set_xticklabels(np.round(np.arange(0 - h_x * dx, Lx, h_x * dx, dtype=np.float64), 1))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(h_y))
        ax1.set_yticklabels(np.round(np.arange(0 - h_y * dx, Ly, h_y * dx, dtype=np.float64), 1))
        plt.xlabel(r'x (m)')
        plt.ylabel(r'y (m)')
        im = plt.imshow(np.real(p_num.T), cmap='viridis')  # , interpolation='none')
        im.set_data(np.real(p_num.T))
        plt.title('It=%3.i, t=%.4f s, FDTD' % (n, n * dt))
        plt.colorbar(format="%0.2f").set_label(label='Pressure (Pa)',fontsize=14)

        ax2 = plt.subplot(212)
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(h_x))
        ax2.set_xticklabels(np.round(np.arange(0 - h_x * dx, Lx, h_x * dx, dtype=np.float64), 1))
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(h_y))
        ax2.set_yticklabels(np.round(np.arange(0 - h_y * dx, Ly, h_y * dx, dtype=np.float64), 1))
        im = plt.imshow(np.real(p_exact.T), cmap='viridis')  # , interpolation='none')
        im.set_data(np.real(p_exact.T))
        plt.xlabel(r'x (m)')
        plt.ylabel(r'y (m)')
        plt.title('It=%3.i, t=%.4f s, Analytic' % (n, n * dt))
        plt.colorbar(format="%0.2f").set_label(label='Pressure (Pa)',fontsize=14)
        plt.pause(10 ** -5)
        plt.tight_layout()
        if n == 10:
            res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                    'case2', 'figures')
            if not os.path.exists(res_path):
                os.makedirs(res_path)
            plt.savefig(os.path.join(res_path, 'modes_nx%i_ny%i.eps' % (nx, ny)),
                        transparent=True, bbox_inches='tight',
                        pad_inches=0)
            plt.savefig(os.path.join(res_path, 'modes_nx%i_ny%i.png' % (nx, ny)),
                        transparent=True, bbox_inches='tight',
                        pad_inches=0)
            plt.savefig(os.path.join(res_path, 'modes_nx%i_ny%i.pdf' % (nx, ny)),
                        transparent=True, bbox_inches='tight',
                        pad_inches=0)


def instatenous_pressure(n, Nt, p, dx, dt, Lx, Ly, case, save_on):
    """
     Plot the instaneous pressure in methods.py if disp_inst_p = True.

    :param n: time iteration index
    :type n: int
    :param Nt: number of iteration time
    :type Nt: int
    :param p: acoustic pressure array (Pa)
    :type p: 2D numpy array of floats.
    :param dx: spatial step (m).
    :type dx: float
    :param dt: time step (s).
    :type dt: float
    :param Lx: continuous length of the domain (in meter) following the x-direction
    :type Lx: float
    :param Ly: continuous length of the domain (in meter) following the y-direction
    :type Ly: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param save_on: boolean to save some snapshots of the wavefront.
    :type save_on: bool
    :return: display (and save snapchots of) the instantaneous pressure in the 2D domain as a function of time.
    """
    # if n > 1 and n % int(Nt / 6.) == 0:
    if n > 1 and n % 20 == 0:
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        plt.imshow(np.real(p.T), cmap='viridis', origin='lower')
        plt.clim(-.5, .5)
        h_x = 40
        h_y = h_x
        ax.xaxis.set_major_locator(ticker.MultipleLocator(h_x))
        ax.set_xticklabels(np.round(np.arange(0 - h_x * dx, Lx, h_x * dx, dtype=np.float64),1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(h_y))
        ax.set_yticklabels(np.round(np.arange(0 - h_y * dx, Ly, h_y * dx, dtype=np.float64),1))

        plt.xlabel(r'x (m)')
        plt.ylabel(r'y (m)')
        plt.title('It=%3.i, t=%.3f s' % (n, n * dt))
        plt.pause(10 ** -50)

        if save_on:
            res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                    'case%i' % case, 'figures')
            if not os.path.exists(res_path):
                os.makedirs(res_path)
            # plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.eps' % (int(dx*100), n)),
            #             transparent=True, bbox_inches='tight', pad_inches=0)
            plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.png' % (int(dx*100), n)),
                        transparent=True, bbox_inches='tight', pad_inches=0)
            # plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.jpg' % (int(dx*100), n)),
            #             transparent=True, bbox_inches='tight', pad_inches=0)
            print 'saving wavefronts'


def instatenous_pressure_fd_method(p, dx, Lx, Ly, case, save_on):
    """
     Plot the instaneous pressure in methods.py if disp_inst_p = True.

    :param n: time iteration index
    :type n: int
    :param Nt: number of iteration time
    :type Nt: int
    :param p: acoustic pressure array (Pa)
    :type p: 2D numpy array of floats.
    :param dx: spatial step (m).
    :type dx: float
    :param dt: time step (s).
    :type dt: float
    :param Lx: continuous length of the domain (in meter) following the x-direction
    :type Lx: float
    :param Ly: continuous length of the domain (in meter) following the y-direction
    :type Ly: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param save_on: boolean to save some snapshots of the wavefront.
    :type save_on: bool
    :return: display (and save snapchots of) the instantaneous pressure in the 2D domain as a function of time.
    """
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    plt.imshow(np.real(p.T), cmap='viridis', origin='lower')
    plt.clim(-10**0, 10**0)
    plt.colorbar()
    # h_x = 40
    # h_y = h_x
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(h_x))
    # ax.set_xticklabels(np.round(np.arange(0 - h_x * dx, Lx, h_x * dx,
    #                                       dtype=np.float64),1))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(h_y))
    # ax.set_yticklabels(np.round(np.arange(0 - h_y * dx, Ly, h_y * dx,
    #                                       dtype=np.float64),1))

    plt.xlabel(r'x (m)')
    plt.ylabel(r'y (m)')
    plt.title('FD method')
    plt.pause(10 ** -50)

    if save_on:
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'figures')
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        # plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.eps' % (int(dx*100), n)),
        #             transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'case0_pressure_1.png' ),
                    transparent=True, bbox_inches='tight', pad_inches=0)

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    fig = plt.figure(2)
    ax = fig.gca(projection='3d')
    X = np.linspace(-Lx / 2., Lx / 2. + dx, p.shape[0])
    Y = np.linspace(-Ly / 2., Ly / 2. + dx, p.shape[1])
    # X = np.arange(-Lx/2., Lx/2. + dx, dx)
    # Y = np.arange(-Ly/2., Ly/2. + dx, dx)
    X, Y = np.meshgrid(X, Y)
    Z = np.real(p.T)
    print np.shape(X), np.shape(Y), np.shape(Z)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.plasma,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-1.01, 1.01)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    if save_on:
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, 'figures')
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        # plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.eps' % (int(dx*100), n)),
        #             transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'case0_pressure_2.png' ),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        # plt.savefig(os.path.join(res_path, 'wf_dx%i_%i.jpg' % (int(dx*100), n)),
        #             transparent=True, bbox_inches='tight', pad_inches=0)
        print 'saving wavefronts'

    plt.show()

def inst_plot(n,p,Nb):
    """

    :param n: time iteration index
    :type n: int
    :param p: acoustic pressure array (Pa)
    :type p: 2D numpy array of floats.
    :param Nb: boundary of the domain (1 if BC, 0 else).
    :type Nb: 2D numpy array of integers.
    :return: show the instantaneous pressure in the 2D domain as a function of time.
    """
    import matplotlib.pyplot as plt
    # instant pressure plot
    if n%10==0:
        plt.figure(1)
        im = plt.imshow(np.real(p.T), cmap='viridis', origin ='lower')
        plt.clim(-.1, .1)
#        im.set_data(np.real(Nb.T))
        im.set_data(np.real(p.T))
        plt.title('It=%i' %(n))
        plt.pause(10**-50)