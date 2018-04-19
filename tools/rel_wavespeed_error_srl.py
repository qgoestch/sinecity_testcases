# -*- coding: utf-8 -*-
##
# \file     rel_wavespeed_error_srl.py
# \title    Calculation of the wave speed error for a Standard Rectilinear Laplacian (Cartesian grid).
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 09 Oct.
##
import numpy as np
import matplotlib.pyplot as plt
import os
base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

def dispersion_relation_k():
    """
         Dispersion relation for 2D Rectilinear grids.
      Calculation of the dispersion relation as a function of the wavenumber.

    :param      c       sound speed, float (m.s - 1).
    :param      dl      spatial step, float (m).
    :param      Ts      time step, float (s).
    :param      theta   azimuthal angle, float (rad).
    :param      k       wavenumber, float (m-1).
    :param      Cn      Courant number, float (m-1).
    :param      kx      wavenumber following the x direction, float (m-1).
    :param      ky      wavenumber followinf the y direction, float (m-1).
    :param      omegSRL_k_theta      numerical wavenumber as a function of k and theta, float (rad.s-1).
    :param      vp      relative phase velocity as a function of k and theta, float.
    :param      omegSRL_k       numerical wavenumber as a function of k, float (rad.s-1).
    :param      vp_k    relative phase velocity as a function of k, float.
    """
    c = 340.0
    dl = 0.16 # worst case for case 1: div. geo.
    Cn = 1. / np.sqrt(2)
    Ts = Cn * dl / c
    kx = np.linspace(-np.pi, np.pi, 5 * 10 ** 2)
    ky = np.linspace(-np.pi, np.pi, 5 * 10 ** 2)
    omegSRL_k = np.zeros((len(kx), len(ky)))

    # Calculation of the relative phase speed error as a function of kx and ky
    vp_k = np.zeros((len(kx), len(ky)))
    for kx_idx, kx_val in enumerate(kx):
        for ky_idx, ky_val in enumerate(ky):
            omegSRL_k[kx_idx, ky_idx] = 2. / Ts * np.arcsin(Cn * np.sqrt(np.sin(dl * kx_val / 2.) ** 2 +
                                                                         np.sin(dl * ky_val / 2.) ** 2))
            vp_k[kx_idx, ky_idx] = omegSRL_k[kx_idx, ky_idx] / (np.sqrt(kx_val ** 2 + ky_val ** 2) * c)

    # Calculation of the relative phase speed error as a function of (kx * dl) and (ky * dl)
    kxdl = np.linspace(-np.pi, np.pi, 5 * 10 ** 2)
    kydl = np.linspace(-np.pi, np.pi, 5 * 10 ** 2)
    omegSRL_kdl = np.zeros((len(kx), len(ky)))
    vp_kdl = np.zeros((len(kx), len(ky)))
    for kxdl_idx, kxdl_val in enumerate(kxdl):
        for kydl_idx, kydl_val in enumerate(kydl):
            omegSRL_kdl[kxdl_idx, kydl_idx] = 2. / Ts * np.arcsin(Cn * np.sqrt(np.sin(kxdl_val / 2.) ** 2 +
                                                                               np.sin(kydl_val / 2.) ** 2))
            vp_kdl[kxdl_idx, kydl_idx] = omegSRL_k[kxdl_idx, kydl_idx] / (np.sqrt((kxdl_val / dl) ** 2 +
                                                                                  (kydl_val / dl) ** 2) * c)
    KX, KY = np.meshgrid(kx, ky)
    fig = plt.figure('Dispersion error (k_x, k_y)')
    ax = fig.add_subplot(111)
    ax.grid(True, which="both", ls=":")
    plt.contourf(KX, KY, np.abs(vp_k).T,
                N=10, norm=None,levels=np.linspace(0.995, 1.00, 50),cmap='viridis',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%0.3f").set_label(label='Rel. phase velocity',fontsize=15)
    CS=plt.contour(KX, KY, np.abs(vp_k).T, levels=np.linspace(0.991, 1.01, 10), colors='k')
    plt.clabel(CS, fontsize=20, inline=3,fmt = '%0.3f',log=True)
    plt.title(r'', fontsize=16)
    plt.xlabel(r'$k_x$ (rad.m$^{-1}$)', fontsize=15)
    plt.ylabel(r'$k_y$ (rad.m$^{-1}$)', fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    fig.tight_layout()
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case1', 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'dispersion.eps'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'dispersion.png'), transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'dispersion.pdf'), transparent=True, bbox_inches='tight', pad_inches=0)

    plt.show()

dispersion_relation_k()

def dispersion_relation_theta_k():
    """
         Dispersion relation for 2D Rectilinear grids.
      Calculation of the dispersion relation as a function of the wavenumber and the angle theta.

    :param      c       sound speed, float (m.s - 1).
    :param      dl      spatial step, float (m).
    :param      Ts      time step, float (s).
    :param      theta   azimuthal angle, float (rad).
    :param      k       wavenumber, float (m-1).
    :param      Cn      Courant number, float (m-1).
    :param      kx      wavenumber following the x direction, float (m-1).
    :param      ky      wavenumber followinf the y direction, float (m-1).
    :param      omegSRL_k_theta      numerical wavenumber as a function of k and theta, float (rad.s-1).
    :param      vp      relative phase velocity as a function of k and theta, float.
    :param      omegSRL_k       numerical wavenumber as a function of k, float (rad.s-1).
    :param      vp_k    relative phase velocity as a function of k, float.
    """
    c = 340.0
    dl = 0.05
    Cn = 0.7 #1. / np.sqrt(2)
    Ts = Cn * dl / c
    theta = np.linspace(0, np.pi, 1 * 10**2)
    k = np.linspace(-np.pi, np.pi, 2 * 10**2)
    omegSRL_k_theta = np.zeros((len(k), len(theta)))
    vp = np.zeros((len(k), len(theta)))

    for k_idx, k_val in enumerate(k):
        for theta_idx, theta_val in enumerate(theta):
            kx = k_val * np.cos(theta_val)
            ky = k_val * np.sin(theta_val)
            omegSRL_k_theta[k_idx, theta_idx] = 2. / Ts * np.arcsin(Cn * np.sqrt(np.sin(dl * kx / 2.) ** 2 +
                                                                                 np.sin(dl * ky / 2.) ** 2))
            vp[k_idx, theta_idx] = omegSRL_k_theta[k_idx, theta_idx] / (np.sqrt(kx ** 2 + ky ** 2) * c)

# dispersion_relation_theta_k()