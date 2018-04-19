# -*- coding: utf-8 -*-
##
# \file     plot_polar_scat.py
# \title    Plot the polar diagramms for each numerical method in case4: scattering.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar), LAUM (Le Mans Université)
# \date     2017, 14 Sep.
##

import numpy as np
import os
from matplotlib import pyplot as plt

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))


def polar_plots(c, radius, f, phi_rcp_grid, p_tlm, p_fdtd, pan,
                f_idx_start, f_idx_end, Mag_Norm_tlm, Mag_Norm_fdtd,
                case, num_method):
    """
      Plot the polar diagramms for a given numerical method.

    :param  c       sound speed, scalar (m.s-1).
    :param  radius  radius of the circular obstacle, float (m).
    :param  f       frequency sequence from the FFT of the numerical simulation, 1d list of floats (Hz).
    :param  phi_rcp_grid    distance of each discrete receivers on the circles on the numerical grid, float (rad).
    :param  p_num   FFT of the time presure signal obtain from the numerical method,
                        3d list [angle,distance,frequency].
    :param  pan     analytic solution for plane waves scattered by a circular obstacle,
                        3d list [angle,distance,frequency].
    :param  f_idx_start     index of the lowest frequency within the range of interest, integer.
    :param  f_idx_end     index of the highest frequency within the range of interest, integer.
    :param  Mag_Norm_tlm  magnitude normalisation of the TLM result in order to fit the analytic solution, float.
    :param  Mag_Norm_fdtd magnitude normalisation of the FDTD result in order to fit the analytic solution, float.
    :param  case    integer that sorts of the saved folders in the results directory.
    :param  num_method      name of the numerical method used in figure title, string.


    """
    f = f.astype(int)
    freq = [f[f_idx_start+5], f[f_idx_start + 10], f[f_idx_start + 20],
            f[f_idx_end - 10], f[f_idx_end - 5], f[f_idx_end - 1]]
    print freq
    f = f.tolist()
    iangle_st = 0
    iangle = -1
    idist_st = int(pan.shape[1] / 2)
    idist = idist_st + 1
    r_idx = 0

    # ==============================================================================
    #   FIGURES: POLAR DIAGRAMS FOR 6 FREQUENCIES
    # ==============================================================================
    plt.figure('Polar patterns: %s' % num_method, figsize=(10, 5.6))
    plt.subplot(231, polar=True)
    i = 0
    j = f.index(freq[i]) - 1
    print 'j=%i' %j
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, f.index(freq[i])]), 'k', lw=1.5)
    print j-f_idx_start+1
    plt.legend(('TLM', 'FDTD', 'Analytic'), bbox_to_anchor=(.5, .3), fontsize=12)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)

    plt.subplot(232, polar=True)
    i = 1
    j = f.index(freq[i]) - 1
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, f.index(freq[i])]), 'k', lw=1.5)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)

    plt.subplot(233, polar=True)
    i = 2
    j = f.index(freq[i]) - 1
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, f.index(freq[i])]), 'k', lw=1.5)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)

    plt.subplot(234, polar=True)
    i = 3
    j = f.index(freq[i]) - 1
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, j-f_idx_start+1]), 'k', lw=1.5)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)

    plt.subplot(235, polar=True)
    i = 4
    j = f.index(freq[i]) - 1
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, j-f_idx_start+1]), 'k', lw=1.5)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)

    plt.subplot(236, polar=True)
    i = 5
    j = f.index(freq[i]) - 1
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_tlm * np.abs(p_tlm[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'rs', markersize=6,
              markeredgewidth=1., markeredgecolor='r', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle:4, idist_st:idist] + np.pi,
              Mag_Norm_fdtd * np.abs(p_fdtd[iangle_st:iangle:4, r_idx, f.index(freq[i])]), 'g+', markersize=10,
              markeredgewidth=1., markeredgecolor='g', markerfacecolor='None')
    plt.polar(phi_rcp_grid[iangle_st:iangle, idist_st:idist] + np.pi,
              np.abs(pan[iangle_st:iangle, r_idx, j-f_idx_start+1]), 'k', lw=1.5)
    plt.title('$f$=%i Hz, $ka=%.2f$' % (freq[i], 2 * np.pi * freq[i] / c * 2 * radius), fontsize=13)
    plt.tight_layout()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'polar_plots.eps'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'polar_plots.png'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'polar_plots.pdf'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.show()
