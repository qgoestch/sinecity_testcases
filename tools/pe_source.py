# -*- coding: utf-8 -*-
##
# \file     pe_source.py
# \title    Parabolic equation source term.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 15 Nov.
##

import numpy as np


def gaussian_source_imp(k, beta, h, h_s, Nx, Ny):
    """
    Spatial Gaussian source that is imposed on the first row of
    the domain. It includes an image source pondered by the ground reflection
    coefficient, see **[blairon_phd2002, Eqs.(2.36)]**.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param beta: **NORMALIZED** admitance used for the boundary condition (kg.s.m-2).
    :type k: float
    :param h: spatial step for both Cartesian directions (m).
    :type h: float
    :param h_s: height of the source (m).
    :type h_s: float
    :param Rp: reflection coefficient of the ground.
    :type Rp: float
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the source pressure to be imposed at the entrance of the domain.
    :rtype: 1D array of floats
    """
    A = 1. / (2. * np.sqrt(2. * np.pi))
    W = np.sqrt(2.) / k     # ouverture de la source
    z_sqr = [(i*h - h_s)**2 for i in range(Ny)]
    z_sqr_img = [(i*h + h_s)**2 for i in range(Ny)]

    # ==============================================================================
    #   Gaussian source with its image source - [chevret_phd1994, Eqs.(4.37), p.59]
    # ==============================================================================
    p_src = np.zeros((Ny + 1, Nx + 1), dtype=np.complex128)
    for i in range(Ny):
        p_src[0, i] = A * np.exp(- z_sqr[i] / W ** 2) + \
                      (1. - beta) / (1. + beta) * \
                      A * np.exp(- z_sqr_img[i] / W ** 2)
    return p_src


def gaussian_source(k, h, h_s, Nx, Ny):
    """
    Spatial Gaussian source that is imposed on the first row of
    the domain, see  - **[chevret_phd1994, Eqs.(4.37), p.59]** or
     **[blairon_phd2002, Eqs.(2.21)-(2.25)]**.

    :param k: wave number: k=2*np.pi*f/c0 (rad.m-1).
    :type k: float
    :param h: spatial step for both Cartesian directions (m).
    :type h: float
    :param h_s: height of the source (m).
    :type h_s: float
    :param Rp: reflection coefficient of the ground.
    :type Rp: float
    :param Nx: length of the domain in number of nodes following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :return: the source pressure to be imposed at the entrance of the domain.
    :rtype: 1D array of floats
    """
    A = 1. / (2. * np.sqrt(2. * np.pi))
    W = np.sqrt(2.) / k     # ouverture de la source
    z_sqr = [(i*h - h_s)**2 for i in range(Ny)]
    z_sqr_img = [(i*h + h_s)**2 for i in range(Ny)]

    # ==============================================================================
    #   Gaussian source with its image source - [chevret_phd1994, Eqs.(4.37), p.59]
    # ==============================================================================
    p_src = np.zeros((Ny + 1, Nx + 1), dtype=np.complex128)
    for i in range(Ny):
        p_src[0, i] = A * np.exp(- z_sqr[i] / W ** 2) + \
                      A * np.exp(- z_sqr_img[i] / W ** 2)
    return p_src
