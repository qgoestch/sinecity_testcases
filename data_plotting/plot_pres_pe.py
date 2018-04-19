# -*- coding: utf-8 -*-
##
# \file     plot_pres_pe.py
# \title    Plots of the absolute pressure given by the parabolic equation.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 16 Nov.
##
import numpy as np
from matplotlib import pyplot as plt


def plot_pressure_pe(p, h_s, freq, dx, dy):
    """
    Display the absolute pressure in the 2D domain of the parabolic equation
    method. The pressure can be pre-processed to show relative or absolute
    levels, or just in Pascal directly from the init_pe_*.

    :param p: pressure from the parabolic equation solver (Pa, dB...).
    :type p: 2D numpy arrays of complexes
    :param h_s: height of the source (m).
    :type h_s: float
    :param freq: frequency of interest (Hz).
    :type freq: float
    :param dx: spatial step for the x directions (m).
    :type dx: float
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :return:
    :rtype:
    """

    z_min = 0.0 * 10**0
    z_max = 140.0 * 10**0

    print np.max(p), np.min(p)

    fig = plt.figure('Parabolic Equation')
    ax = fig.add_subplot(111)
    # plt.grid()
    x5 = [i * dx for i in range(p.shape[0])]
    y5 = [i * dy for i in range(p.shape[1])]
    X5, Y5 = np.meshgrid(x5, y5)
    plt.contourf(X5, Y5, np.abs(p[:, :].T),
                N=40, norm=None, levels=np.linspace(z_min, z_max, 150),
                cmap='viridis', interpolation='nearest', origin='lower')
    plt.colorbar(format="%.i").\
        set_label(label='$20log_{10}(|p(x,y)|/2e^{-5})$ (dB-SPL)',
                  fontsize=12)
    plt.title(r'$h_s$=%.2fm, f=%.2f Hz' % (h_s, freq), fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('y (m)', fontsize=12)
    plt.axis('equal')
    fig.tight_layout()

    plt.show()