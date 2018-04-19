# -*- coding: utf-8 -*-
##
# \file     plot_pres_slope_rot.py
# \title    Plots of the absolute pressure given by the analytic solution (Rasmussen).
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 03 Oct.
##
import numpy as np
from matplotlib import pyplot as plt
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def plot_pres_slope_with_rot(l1, l2, h, res03, res13, h_s,freq, alpha, dx, model_name):
    """
    Two figures that gives the absolute pressure obtained from the analytic
    solution (Rasmussen).

    :param      l1      discretized length of the horizontal part, array.
    :param      l2      discretized length of the upward-sloping part, array.
    :param      h       discretized heights of both parts, array.
    :param      res03   absolute value of the total pressure: incident+diffracted for the flat part, float (Pa).
    :param      res13   absolute value of the total pressure: incident+diffracted after corner, float (Pa).
    :param      res12   absolute value of the diffracted pressure after corner, float (Pa).
    :param      freq    frequency of the source signal, float (Hz).
    :param      alpha   angle of the upward-sloping part, float (rad).
    :param      model_name      name of the analytic mmodel, string.

    :param  z_min   minimal magnitude for the scale of the colorbar, float.
    :param  z_max   maximal magnitude for the scale of the colorbar, float.
    """
    z_min=1.0 * 10**-5
    z_max=0.055
    l3 = [i + l1[-1] for i in l2]
    l = np.concatenate((l1, l3), axis=0)
    res = np.concatenate((res03,res13), axis=0)
    from scipy import ndimage
    rotate_res13 = ndimage.rotate(res13.astype(float), alpha)

    l4 = [i * np.cos(alpha * np.pi / 180.) + l1[-1] for i in l2]
    l44 = np.concatenate((l1, l4), axis=0)

    new_shape = (res03.shape[0], rotate_res13.shape[1])  # This will be some pre-determined size
    res03_reshape = np.zeros(new_shape)
    print np.shape(res03), np.shape(rotate_res13)
    res03_reshape[:, :res03.shape[1]] = res03
    print np.shape(res03), np.shape(res03_reshape), np.shape(rotate_res13)
    rotate_res = np.concatenate((res03_reshape, rotate_res13), axis=0)

    dist_del_idx = int(round((h[-1] * np.sin(alpha * np.pi / 180.)) / dx))
    start_del_idx = int(round(l1[-1] / dx)) - 1
    end_del_idx = start_del_idx + dist_del_idx + 1
    rotate_res_del = np.delete(rotate_res, range(start_del_idx, end_del_idx), 0)

    fig = plt.figure('Analytic model, ROT', figsize=(6., 3.))
    ax = fig.add_subplot(111)
    # plt.grid()
    x5 = [i * dx for i in range(rotate_res_del.shape[0])]
    y5 = [i * dx for i in range(rotate_res_del.shape[1])]
    X5, Y5 = np.meshgrid(x5, y5)
    plt.contourf(X5, Y5, rotate_res_del[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'$h_s$=%.2fm, f=%.2f Hz, $\alpha$ = %.0f$^{\circ}$' %(h_s, freq, alpha),
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
    plt.savefig(os.path.join(res_path, 'analytic_%iHz_%ideg.eps' % (int(freq), int(alpha))),
                            transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'analytic_%iHz_%ideg.png' % (int(freq), int(alpha))),
                            transparent=True, bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'analytic_%iHz_%ideg.pdf' % (int(freq), int(alpha))),
                            transparent=True, bbox_inches='tight', pad_inches=0)
    plt.show()
