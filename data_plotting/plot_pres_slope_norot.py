# -*- coding: utf-8 -*-
##
# \file     plot_pres_slope_norot.py
# \title    Plots of the absolute pressure given by numerical models.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar), LAUM (Le Mans Université)
# \date     2017, 05 Oct.
##
import numpy as np
from matplotlib import pyplot as plt
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def plot_pres_slope_without_rotation(res03, res13, alpha, h_s, freq, dx, n_idx_step, num_meth):
    """
    Two figures that gives the absolute pressure obtained from the analytic
    solution (Rasmussen).

    :param res03: absolute value of the total pressure: incident+diffracted for the flat part (Pa).
    :type res03: array of floats
    :param res13: absolute value of the total pressure: incident+diffracted after corner (Pa).
    :type res13: array of floats
    :param alpha: angle of the upward-sloping part (rad).
    :type alpha: float
    :param h_s: source heights (m).
    :type h_s: float
    :param freq: frequency of the source signal (Hz).
    :type freq: float
    :param dx: spatial step (m).
    :type dx: float
    :param n_idx_step: number of indexes skipped during the FFT calculation
    :type n_idx_step: int
    :param num_meth: name of the numerical method
    :type num_meth: str

    :param z_min: minimal magnitude for the scale of the colorbar, float.
    :param z_max: maximal magnitude for the scale of the colorbar, float.

    :return: 2D plot of the acoustic pressure for a given frequency
    """

    print np.shape(res03), np.shape(res13)
    res = np.concatenate((res03, res13))
    z_min = 1.0 * 10**-10
    z_max = 5.0 * 10**-2 #2.2 * np.average(np.abs(res))

    # print np.abs(res[int(res.shape[0]/2):res.shape[0], :].T)

    fig = plt.figure('Numerical model, NO ROT, %s' % num_meth, figsize=(6., 3.))
    ax = fig.add_subplot(111)
    # plt.grid()
    x5 = [i * n_idx_step * dx for i in range(res.shape[0])]
    y5 = [i * n_idx_step * dx for i in range(res.shape[1])]
    X5, Y5 = np.meshgrid(x5, y5)
    plt.contourf(X5, Y5, np.abs(res[:, :].T),
                N=40, norm=None, levels=np.linspace(z_min, z_max, 150), cmap='viridis',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2e").set_label(label='$|p(x,y)|$ (Pa)', fontsize=12)
    plt.title(r'$h_s$=%.2fm, f=%.2f Hz, $\alpha$ = %.0f$^{\circ}$' % (h_s, freq, alpha),
              fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)', fontsize=12)
    # plt.axis('equal')
    fig.tight_layout()
    plt.ylim(0, 50)
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results', 'case5', 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, '%s_%iHz_%ideg.eps' % (num_meth, int(freq), int(alpha))),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, '%s_%iHz_%ideg.png' % (num_meth, int(freq), int(alpha))),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, '%s_%iHz_%ideg.pdf' % (num_meth, int(freq), int(alpha))),
                transparent=True, bbox_inches='tight', pad_inches=0)
    plt.show()
