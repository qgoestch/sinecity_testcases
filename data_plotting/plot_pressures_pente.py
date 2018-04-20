# -*- coding: utf-8 -*-
##
# \file     plot_pressures_pente.py
# \title    Plots of the absolute pressure given by the analytic solution (Rasmussen).
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 29 May
##
import numpy as np
from matplotlib import pyplot as plt
def plot_pressure_heights(l1,l2,h,res03,res13,h_s,freq,alpha,model_name):
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
    z_min=0.00
    z_max=0.25
    X1, Y1 = np.meshgrid(l1[:],h)
    X2, Y2 = np.meshgrid(l2[:],h)
    print np.shape(X1),np.shape(Y1),np.shape(res03[:,:].T)
    print np.shape(X2), np.shape(Y2),np.shape(res13[:,:].T)
    # fig = plt.figure('%s 1' %model_name)
    # ax = fig.add_subplot(311)
    # plt.grid()
    # plt.contourf(X1,Y1, np.transpose(res03[:,:]),\
    #             N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',\
    #             interpolation='nearest', origin='lower')
    # plt.colorbar(format="%1.2e").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.ylabel('y (m)',fontsize=12)
    # fig.tight_layout()
    # ax = fig.add_subplot(312)
    # plt.grid()
    # plt.contourf(X2,Y2, np.transpose(res13[:,:]),\
    #             N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',\
    #             interpolation='nearest', origin='lower')
    # plt.colorbar(format="%1.2e").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.ylabel('y (m)',fontsize=12)
    # fig.tight_layout()
    # ax = fig.add_subplot(313)
    # plt.grid()
    # plt.contourf(X2,Y2, np.transpose(res12[:,:]),\
    #             N=40, norm=None,levels=np.linspace(z_min, 0.008,150),cmap='magma_r',\
    #             interpolation='nearest', origin='lower')
    # plt.colorbar(format="%1.2e").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.xlabel('x (m)',fontsize=12)
    # plt.ylabel('y (m)',fontsize=12)
    # plt.yticks(fontsize=15)
    # fig.tight_layout()

    z_min=0.00
    z_max=0.50

    fig = plt.figure('%s 2' %model_name)
    ax = fig.add_subplot(211)
    plt.grid()
    plt.contourf(X1,Y1, res03[:,:].T,\
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',\
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'Flat part: $h_s$=%.2fm, f=%.2f Hz' %(h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    fig.tight_layout()
    X, Y = np.meshgrid(l2[:-1],h)
    ax = fig.add_subplot(212)
    plt.grid()
    plt.contourf(X2,Y2, res13[:,:].T,\
                N=40, norm=None,levels=np.linspace(z_min, z_max,150),cmap='viridis_r',\
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.2f").set_label(label='$|p(x,y)|$ (Pa)',fontsize=12)
    plt.title(r'After corner: $\alpha=%i^\circ$, f=%.2f Hz' %(alpha,freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    fig.tight_layout()
    plt.show()
