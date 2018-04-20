# -*- coding: utf-8 -*-
##
# \file     plot_pressures_pente.py
# \title    Plots of the absolute pressure given by the analytic solution (Rasmussen).
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 03 Oct.
##
import numpy as np
from matplotlib import pyplot as plt
def plot_pressure_slope(l1,l2,h,res03,res13,h_s,freq,alpha, dx, model_name):
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
    """
    X1, Y1 = np.meshgrid(l1[:],h)
    X2, Y2 = np.meshgrid(l2[:],h)
    z_min=0.0001
    z_max=0.10
    fig = plt.figure('Analytic model, %s 2' %model_name)
    ax = fig.add_subplot(311)
    plt.grid()
    plt.contourf(X1,Y1, res03[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'Flat part: $h_s$=%.2fm, f=%.2f Hz' %(h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    plt.axis([0, l1[-1], 0, h[-1]])
    plt.xlim(0, l1[-1])
    fig.tight_layout()
    ax = fig.add_subplot(312)
    plt.grid()
    plt.contourf(X2,Y2, res13[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'After corner: $\alpha=%i^\circ$, f=%.2f Hz' %(alpha,freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    fig.tight_layout()

    l3 = [i + l1[-1] for i in l2]
    l = np.concatenate((l1, l3), axis=0)
    res = np.concatenate((res03,res13), axis=0)
    from scipy import ndimage
    rotate_res13 = ndimage.rotate(res13.astype(float), alpha)
    print np.shape(rotate_res13), np.shape(res13), rotate_res13.shape[0], rotate_res13.shape[1]
    # print np.shape(res), len(X), len(Y),l[-1]

    X3, Y3 = np.meshgrid(range(rotate_res13.shape[0]), range(rotate_res13.shape[1]))
    ax = fig.add_subplot(313)
    plt.grid()
    plt.contourf(X3,Y3, rotate_res13[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'After corner: $\alpha=%i^\circ$, f=%.2f Hz' %(alpha,freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    fig.tight_layout()

    X, Y = np.meshgrid(range(res.shape[0]), range(res.shape[1]))
    fig = plt.figure('Analytic model, %s 3' %model_name)
    ax = fig.add_subplot(111)
    plt.grid()
    plt.contourf(X,Y, res[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'Flat part: $h_s$=%.2fm, f=%.2f Hz' %(h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    plt.xlim(0, l3[-1])
    fig.tight_layout()

    l4 = [i * np.cos(alpha * np.pi / 180.) + l1[-1] for i in l2]
    l44 = np.concatenate((l1, l4), axis=0)

    new_shape = (res03.shape[0], rotate_res13.shape[1])  # This will be some pre-determined size
    res03_reshape = np.zeros(new_shape)
    print np.shape(res03), np.shape(rotate_res13)
    res03_reshape[:, :res03.shape[1]] = res03
    print np.shape(res03), np.shape(res03_reshape), np.shape(rotate_res13)
    rotate_res = np.concatenate((res03_reshape, rotate_res13), axis=0)
    fig = plt.figure('Analytic model, %s 4' %model_name)
    ax = fig.add_subplot(111)
    plt.grid()
    X4, Y4 = np.meshgrid(range(rotate_res.shape[0]), range(rotate_res.shape[1]))
    plt.contourf(X4,Y4, rotate_res[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'Flat part: $h_s$=%.2fm, f=%.2f Hz' %(h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    # plt.xlim(0, l44[-1])
    fig.tight_layout()

    dist_del_idx = int(round((h[-1] * np.sin(alpha * np.pi / 180.)) / dx))
    start_del_idx = int(round(l1[-1] / dx)) - 1
    end_del_idx = start_del_idx + dist_del_idx + 1
    rotate_res_del = np.delete(rotate_res, range(start_del_idx, end_del_idx), 0)

    fig = plt.figure('Analytic model, %s 5' %model_name)
    ax = fig.add_subplot(111)
    plt.grid()
    X5, Y5 = np.meshgrid(range(rotate_res_del.shape[0]), range(rotate_res_del.shape[1]))
    plt.contourf(X5, Y5, rotate_res_del[:,:].T,
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'Flat part: $h_s$=%.2fm, f=%.2f Hz' %(h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    plt.axis('equal')
    # plt.xlim(0, l44[-1])
    fig.tight_layout()

    plt.show()