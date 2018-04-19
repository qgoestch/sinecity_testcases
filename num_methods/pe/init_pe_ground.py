# -*- coding: utf-8 -*-
##
# \file     init_pe_ground.py
# \title    Definition of the numerical parameters for the PE method.
#           The updated scheme that gives the pressure
#           at each time iteration is defined in the upd_pe.py module.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 14 Nov.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                    os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

print 'current work dir = %s' % os.getcwd()
print 'base path        = %s' % base_path
pe_core_path = os.path.join(base_path, 'pe_core')
site.addsitedir(pe_core_path)
from upd_pe import pe_scheme_pade_22_full, \
                    pe_pade_one_one, pe_pade_two_two, \
                    pe_solver

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from miki_imp_model import Miki
from pe_source import gaussian_source, gaussian_source_imp
from pe_absorbing_layer import abs_lay_top, abs_lay_top_1, abs_lay_bottom_top

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_pres_pe import plot_pressure_pe


def pe_init_impgr(delta_x, d_sr, h_s, h_r, freq, rho, c, sigma, case,
                  free_field, disp):
    """
    Definition of the numerical domain for the application of the parabolic
    equation to the ground reflection case, with an impedance ground defined
    using Miki's impedance model.

    :param delta_x: spatial step for the x directions (m).
    :type delta_x: float
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param freq: frequency of interest (Hz).
    :type freq: float
    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param sigma: pecific airflow resistivity (kNm-4s==CGS).
    :type sigma: float
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param free_field: the domain is enlarged
    :type free_field: bool
    :param disp: display the pressure inside the numerical domain.
    :type disp: bool

    :return: plot and/or save the pressure arrays for post processing.
    """
    # =========================================================================
    #   Grid parameters
    # =========================================================================
    k = 2. * np.pi * freq / c       # wave number (rad.m-1)
    Lx = 1.02 * d_sr[-1]            # length of the horizontal direction (m)
    if not free_field:
        Ly = 10. * max(h_s, h_r)    # length of the vertical direction (m)
    else:
        Ly = 2. * 10. * max(h_s, h_r)
    Nx = np.int(np.round(Lx / delta_x))     # number of node following x axis
    Ny = Nx                         # number of node following y axis
    x = np.linspace(0, Lx, Nx + 1)  # discretized x axis
    y = np.linspace(0, Ly, Ny + 1)  # discretized y axis
    dx = np.float64(x[1] - x[0])    # spatial step for the x direction
    dx = round(dx, 5)
    dy = np.float64(y[1] - y[0])    # spatial step for the y direction
    dy = round(dy, 5)

    print '==============================================================='
    print '               PE above an impedance ground                    '
    print '---------------------------------------------------------------'
    if free_field:
        print 'Free field calculation.'
    print 'Frequency:    f=%0.0f Hz' % freq
    print 'Resistivity:  sigma=%i' % sigma
    print 'DIMENSIONS:   Nx=%i cells ; Ny=%i cells; Lx=%g m; Ly=%g m.' \
                            % (Nx, Ny, Lx, Ly)
    print 'SPATIAL-STEP: dx=%g m, dy=%g m.' % (dx, dy)
    print 'Source height: h_s = %.1fm.' % h_s
    print 'Receivers: x_ref = %0.1fm ; x_meas = %0.1fm. \n' % (d_sr[0], d_sr[1])

    # =========================================================================
    #   Pressure variables
    # =========================================================================
    p_i1j = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)  # final pres. array
    p_saved = np.zeros((len(d_sr)), dtype=np.complex128)

    # =========================================================================
    #   Location of the source and receiver(s) in number of nodes
    # =========================================================================
    x_rcv = [int(round(ii / dx)) for ii in d_sr]
    if free_field:
        abs_lay_len = abs(Ly / 2. - 0.6)
        sigma = 1.0 * 10 ** - 6
        y_rcv = int(round(Ly / 2. / dx))
        h_s = Ly / 2.
    else:
        abs_lay_len = abs(Ly - 1.1 * max(h_r, h_s))
        y_rcv = int(round(h_r / dx))

    # =========================================================================
    #   Surface impedance for the ground boundary condition
    # =========================================================================
    z_miki, k_miki = Miki(-1, freq, sigma, rho, c)
    b_miki = rho * c / z_miki

    # =========================================================================
    #   Matrices initialization for homogeneous propagation
    # =========================================================================
    # Aij, Bij = pe_pade_two_two(k, b_miki, dx, dy, Nx, Ny)
    Aij, Bij = pe_pade_one_one(k, b_miki, dx, dy, Nx, Ny)

    # =========================================================================
    #   Source initialization
    # =========================================================================
    p_ij = gaussian_source(k, dy, h_s, Nx, Ny)
    # p_ij = gaussian_source_imp(k, b_miki, dy, h_s, Nx, Ny)

    # =========================================================================
    #   Calculation of the pressure
    # =========================================================================
    norm_distance = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    for i in range(1, Nx):
        for j in range(1, Ny):
            norm_distance[i, j] = 1. / np.sqrt((i * dx) ** 2 +
                                               (h_s - j * dy) ** 2)

    for i in range(0, Nx):
        # print Nx - i
        # p_ij[i + 1, :] = pe_scheme_pade_22_full(k, b_miki, dx, dy,
        #                                          p_ij[i, :], Nx, Ny)

        # Matrices multiplication and inversion
        p_ij[i + 1, :] = pe_solver(Aij, Bij, p_ij[i, :])

        # Definition of the absorbing layer on top of the domain
        if free_field:
            p_ij[i + 1, :] = abs_lay_bottom_top(p_ij[i + 1, :], dy, Ny,
                                                Ly - abs_lay_len)
        else:
            p_ij[i + 1, :] = abs_lay_top(p_ij[i + 1, :], dy, Ny,
                                         Ly - abs_lay_len)

        # Normalization of the pressure as a function of the distance
        p_i1j[i + 1, :] += np.abs(p_ij[i + 1, :]) * norm_distance[i + 1, :]

    # print np.max(np.abs(p_i1j)), np.min(np.abs(p_i1j))

    # =========================================================================
    #   Display the pressure inside the numerical domain
    # =========================================================================
    if disp:
        plot_pressure_pe(10. * np.log10(np.abs(p_i1j)**2 / (2. * 10 ** -5)**2),
                         h_s, freq, dx, dy)

    # =========================================================================
    #   Save the results in .npy or .npz (compressed) formats
    # =========================================================================
    for d in range(len(d_sr)):
        p_saved[d] = p_i1j[x_rcv[d], y_rcv]

    import os
    res_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                            'results', 'case%i' % case, 'pe')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    if free_field:
        np.save(os.path.join(res_path, 'p_%s_%iHz.npy' % ('f', int(freq))),
                p_saved)
    else:
        np.save(os.path.join(res_path, 'p_%s_%iHz_%icgs.npy'
                             % ('t', int(freq), int(sigma))), p_saved)
