# -*- coding: utf-8 -*-
##
# \file     plot_decreases_geospr.py
# \title    Theoretical and geometrical decreases.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 12 Oct.
##
import os
import matplotlib.pyplot as plt
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )

def plot_decreases(theo_axial_fdtd, pmax_axial_fdtd, theo_diag_fdtd, pmax_diag_fdtd,
                   theo_axial_tlm, pmax_axial_tlm, theo_diag_tlm, pmax_diag_tlm, d_sr):
    """
    :param theo_axial_fdtd: theoretical pressure decrease normalized at 1m distance, float (Pa).
    :param pmax_axial_fdtd: pressure decrease normalized at 1m distance from FDTD simaulations
                            on axial nodes, float (Pa).
    :param theo_diag_fdtd: theoretical pressure decrease normalized at 1m distance, float (Pa).
    :param pmax_diag_fdtd: pressure decrease normalized at 1m distance from FDTD simulations
                           on diagonal nodes, float (Pa).
    :param theo_axial_tlm: theoretical pressure decrease normalized at 1m distance, float (Pa).
    :param pmax_axial_tlm: pressure decrease normalized at 1m distance from TLM simulations
                           on axial nodes, float (Pa).
    :param theo_diag_tlm: theoretical pressure decrease normalized at 1m distance, float (Pa).
    :param pmax_diag_tlm: pressure decrease normalized at 1m distance from TLM simulations
                           on diagonal nodes, float (Pa).
    :param d_sr: distance the aquisition points and the source, float (m).
    :return: 4 graphs that show the comparison between theoretical and numerical pressure decrease for each method and
             each spatial step
    """

    plt.figure('Amplitude decrease axial FDTD')
    plt.plot(d_sr[1:], theo_axial_fdtd[0, :], 'k', lw=2)
    plt.plot(d_sr[1:], pmax_axial_fdtd[0, :], 'ro')
    plt.plot(d_sr[1:], pmax_axial_fdtd[1, :], 'b*')
    plt.plot(d_sr[1:], pmax_axial_fdtd[2, :], 'ys')
    plt.plot(d_sr[1:], pmax_axial_fdtd[3, :], 'm+')
    plt.plot(d_sr[1:], pmax_axial_fdtd[4, :], 'gx')
    plt.grid()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Distance (m)', fontsize=16)
    plt.ylabel('Nomalized pressure max.', fontsize=14)
    plt.title('FDTD axial')
    plt.legend(('Theory', 'FDTD, h1', 'h2', 'h3', 'h4', 'h5'), loc=1, fontsize=12)

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case1', 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'pres_dec_fdtd.eps'), transparent=True, bbox_inches='tight',
                pad_inches=0)
    plt.savefig(os.path.join(res_path, 'pres_dec_fdtd.png'), transparent=True, bbox_inches='tight',
                pad_inches=0)
    plt.savefig(os.path.join(res_path, 'pres_dec_fdtd.pdf'), transparent=True, bbox_inches='tight',
                pad_inches=0)
    plt.show()

    plt.figure('Amplitude decrease diagonal FDTD')
    plt.plot(d_sr[1:], theo_diag_fdtd[0, :], 'k', lw=2)
    plt.plot(d_sr[1:], pmax_diag_fdtd[0, :], 'ro')
    plt.plot(d_sr[1:], pmax_diag_fdtd[1, :], 'b*')
    plt.plot(d_sr[1:], pmax_diag_fdtd[2, :], 'ys')
    plt.plot(d_sr[1:], pmax_diag_fdtd[3, :], 'm+')
    plt.plot(d_sr[1:], pmax_diag_fdtd[4, :], 'gx')
    plt.grid()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Distance (m)', fontsize=16)
    plt.ylabel('Nomalized pressure max.', fontsize=14)
    plt.title('FDTD diagonal')
    plt.legend(('Theory', 'FDTD, h1', 'h2', 'h3', 'h4', 'h5'), loc=1, fontsize=12)

    plt.figure('Amplitude decrease axial TLM')
    plt.plot(d_sr[1:], theo_axial_tlm[0, :], 'k', lw=2)
    plt.plot(d_sr[1:], pmax_axial_tlm[0, :], 'ro')
    plt.plot(d_sr[1:], pmax_axial_tlm[1, :], 'b*')
    plt.plot(d_sr[1:], pmax_axial_tlm[2, :], 'ys')
    plt.plot(d_sr[1:], pmax_axial_tlm[3, :], 'm+')
    plt.plot(d_sr[1:], pmax_axial_tlm[4, :], 'gx')
    plt.grid()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Distance (m)', fontsize=16)
    plt.ylabel('Nomalized pressure max.', fontsize=14)
    plt.title('TLM axial')
    plt.legend(('Theory', 'TLM, h1', 'h2', 'h3', 'h4', 'h5'), loc=1, fontsize=12)

    plt.figure('Amplitude decrease diagonal TLM')
    plt.plot(d_sr[1:], theo_diag_tlm[0, :], 'k', lw=2)
    plt.plot(d_sr[1:], pmax_diag_tlm[0, :], 'ro')
    plt.plot(d_sr[1:], pmax_diag_tlm[1, :], 'b*')
    plt.plot(d_sr[1:], pmax_diag_tlm[2, :], 'ys')
    plt.plot(d_sr[1:], pmax_diag_tlm[3, :], 'm+')
    plt.plot(d_sr[1:], pmax_diag_tlm[4, :], 'gx')
    plt.grid()
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel('Distance (m)', fontsize=16)
    plt.ylabel('Nomalized pressure max.', fontsize=14)
    plt.title('TLM diagonal')
    plt.legend(('Theory', 'TLM, h1', 'h2', 'h3', 'h4', 'h5'), loc=1, fontsize=12)
    plt.show()