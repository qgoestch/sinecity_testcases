# -*- coding: utf-8 -*-
##
# \file     analytic_up_slop.py
# \title    Rasmussen's analytic solution for the diffraction from a corner.
#           See [rasmussen_jsv1985] and [blairon_phd2002, Sec.1.2.1, p.15; Sec.3.1 p.92-99].
# \author   Pierre Chobeau (rewritten from N. Blairon - ECL 2002)
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 26 Sep.
##
import numpy as np
import scipy.special as sp
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))
data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
print data_plotting_path
from plot_pres_slope_rot import plot_pres_slope_with_rot


def main_analytic_slope(freq, h_s, l1max, l2max, hmax, alpha, dx):
    """
      Calculation of the analytic pressure for each point of both
            discretized parts, i.e. horizontal and upward-slopping,
            based on [rasmussen_jsv1985].

    :param      freq    frequency of the source signal, float (Hz).
    :param      h_s     source heights, float (m).
    :param      l1_max  length of the horizontal part, float (m).
    :param      l2_max  length of the upward-sloping part, float (m).
    :param      hmax    height of the the two parts, float (m).
    :param      alpha   angle of the upward-sloping part, float (rad).
    :param      dx      spatial step, float (m).

    :param  c       sound speed (m.s-1), scalar.
    :param  lamb    wave length (m), scalar.
    :param  k       wave number (rad.m-1), scalar.
    :param  a       calculation parameter used in the v function, scalar.
    :param  alpha   angle of the upward-sloping part (rad), scalar.
    :param  gamma   angle between the two parts following a counterclockwise
                        direction (rad), scalar.
    :param  nu      angular variable used in the v function (rad), scalar.

    :param  dim1    number of points for the horizontal part, scalar.
    :param  dim2    number of points for the upward-sloping part, scalar.
    :param  dimh    number of points for the height of both parts, scalar.

    :param  dx1     spatial step follwing the length of the horizontal part, scalar.
    :param  dx2     spatial step follwing the length of the upward-sloping part, scalar.
    :param  dxh     spatial step follwing the height of both parts, scalar.

    :param  l1      discretized length of the horizontal part, array.
    :param  l2      discretized length of the upward-sloping part, array.
    :param  h       discretized heights of both parts, array.

    :param  r1,r2,r3,r4 the four paths, scalar.
    :param  p0      pressure field of the horizontal part, complex.
    :param  pdiff0  diffracted pressure field of the horizontal part, complex.
    :param  p1      pressure field of the upward-slopping part, complex.
    :param  pdiff1  diffracted pressure field of the upward-slopping part, complex.
    :param  res01   absolute value of the incident pressure for the flat part, float (Pa).
    :param  res02   absolute value of the diffracted pressure for the flat part, float (Pa).
    :param  res03   absolute value of the total pressure: incident+diffracted for the flat part, float (Pa).
    :param  res11   absolute value of the incident pressure for the upward-slopping part, float (Pa).
    :param  res12   absolute value of the diffracted pressure for the upward-slopping part, float (Pa).
    :param  res13   absolute value of the total pressure: incident+diffracted for the upward-slopping part, float (Pa).
    :return The absolute values of the pressure fields for each point of
            each domain in the 2D-arrays res01 to res13.
    """
    c = 340.
    lamb = c / freq
    k = 2. * np.pi / lamb
    alpha = np.pi * alpha / 180.
    gamma = np.arctan(h_s / l1max)
    nu = 2. - (np.pi + alpha) / np.pi

    l1 = np.arange(0, l1max + dx, dx)
    l2 = np.arange(0, l2max + dx, dx)
    h = np.arange(dx, hmax + dx, dx)
    # ==============================================================================
    #    # ------- Horizontal part (dim1*dimh), before corner ------- #
    # ==============================================================================
    res01 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    res02 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    res03 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    for mm in range(1, len(l1)):
        #        print dim1-mm
        beta0 = np.arctan(h / (l1max - l1[mm]))
        beta = np.pi - alpha - beta0

        r1 = np.sqrt(l1[mm] ** 2 + (h - h_s) ** 2)
        r2 = np.sqrt(l1[mm] ** 2 + (h + h_s) ** 2)

        p0 = (1. / 4) * np.sqrt(2. / (np.pi * k * r1)) * (np.exp(1j * k * r1)) \
             + (1. / 4) * np.sqrt(2. / (np.pi * k * r2)) * (np.exp(1j * k * r2))

        # ------- Diffracted part ------- #
        r0diff = np.sqrt(h_s ** 2 + l1max ** 2)
        r1diff = np.sqrt((l1max - l1[mm]) ** 2 + h ** 2)
        gr1diff = r0diff + r1diff
        pdiff0 = (1. / 4) * np.sqrt(2. / (np.pi * k * gr1diff)) * np.exp(1j * k * (gr1diff)) \
                 * (v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta + gamma, nu, len(h), k)
                    + v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta - gamma, nu, len(h), k))

        res01[mm - 1, :] = np.abs(p0)
        res02[mm - 1, :] = np.abs(pdiff0)
        res03[mm - 1, :] = np.abs(p0 + pdiff0)

    # ==============================================================================
    #      # ------- Upward-sloping part (dim2*dimh) after corner ------- #
    # ==============================================================================
    res11 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    res12 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    res13 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    for mm in range(1, len(l2)):
        #        print dim2-mm
        beta = np.asarray([np.pi / 2 for i in range(len(h))])
        if mm > 2:
            beta = np.arctan(h / l2[mm])

        r0diff = np.sqrt(h_s ** 2 + l1max ** 2)
        r1diff = np.sqrt(l2[mm] ** 2 + h ** 2)
        gr1diff = r0diff + r1diff
        r1 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha + gamma + beta))
        r2 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha - gamma + beta))
        r3 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha + gamma - beta))
        r4 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha - gamma - beta))
        h4 = l2[mm] * np.tan(alpha - gamma)
        h3 = l2[mm] * np.tan(alpha + gamma)
        p1 = np.zeros((len(h)), dtype=np.complex128)
        for ll in range(len(h)):
            p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r1[ll])) * (np.exp(1j * k * r1[ll])) \
                      + (1. / 4) * np.sqrt(2. / (np.pi * k * r2[ll])) * (np.exp(1j * k * r2[ll]))
            if h[ll] <= h3:
                p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r3[ll])) * (np.exp(1j * k * r3[ll]))
            if h[ll] <= h4:
                p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r4[ll])) * (np.exp(1j * k * r4[ll]))

        # ------- Diffracted part ------- #
        pdiff1 = (1. / 4) * np.sqrt(2. / (np.pi * k * gr1diff)) * np.exp(1j * k * (gr1diff)) \
                 * (v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta + gamma, nu, len(h), k)
                    + v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta - gamma, nu, len(h), k))

        res11[mm - 1, :] = np.abs(p1)
        res12[mm - 1, :] = np.abs(pdiff1)
        res13[mm - 1, :] = np.abs(p1 + pdiff1)

    model_name = 'Rasmussen upward slope'
    plot_pres_slope_with_rot(l1, l2, h, res03, res13, h_s, freq, alpha * 180. / np.pi, dx, model_name)
    return res03, res13


def fp(x):
    #   real, dimension(dimh) :: x
    #   real, dimension(dimh) :: rx
    rx = np.sqrt(x * 2 / np.pi)
    s_fresnel, c_fresnel = sp.fresnel(rx)
    return - 2 * 1j * np.sqrt(x) * np.exp(-1j * x) * np.sqrt(np.pi / 2.) * \
            (.5 - c_fresnel + 1j * (.5 - s_fresnel))


def v(a, b, th, nu, dimh, k):
    #    real, b
    #    real, dimension(dimh) :: a, th
    #    real, dimension(dimh) :: vp, vm, v
    #    integer, dimension(dimh) :: nplus, nmoins
    vp = np.zeros((dimh), dtype=np.complex128)
    vm = np.zeros((dimh), dtype=np.complex128)
    nmoins = np.zeros((dimh), dtype=np.float64)
    nplus = np.ones((dimh), dtype=np.float64)
    for th_idx, th_val in enumerate(th):
        if (th_val - nu * np.pi + np.pi) <= 0:
            nplus[th_idx] = 0
        nmoins[th_idx] = 1
        if (th_val - nu * np.pi - np.pi) <= 0:
            nmoins[th_idx] = 0
        if (th_val + nu * np.pi - np.pi) <= 0:
            nmoins[th_idx] = -1
        vp[th_idx] = -  (np.exp(1j * np.pi / 4)) * (2 * np.pi * k * a[th_idx] * b) ** (-.5) / (2 * nu) * \
                     (np.tan((np.pi + th_val) / (2 * nu))) ** (-1) * \
                     fp(2 * k * a[th_idx] * (np.cos((2 * nplus[th_idx] * nu * np.pi - th_val) / 2)) ** 2)
        vm[th_idx] = -  (np.exp(1j * np.pi / 4)) * (2 * np.pi * k * a[th_idx] * b) ** (-.5) / (2 * nu) * \
                     (np.tan((np.pi - th_val) / (2 * nu))) ** (-1) * \
                     fp(2 * k * a[th_idx] * (np.cos((2 * nmoins[th_idx] * nu * np.pi - th_val) / 2)) ** 2)
    return vp + vm


if __name__ == '__main__':
    """  
    Test the analytic solution such as in **[blairon_phd2002, Sec.3.3.1, p.87]**
    :param freq: frequency of the source signal (Hz), scalar.
    :param h_s: height of the source (m) located at the entrance of the
                        flat part, scalar.
    :param l1max: length of the horizontal part (m), scalar.
    :param l2max: length of the upward-sloping part (m), scalar.
    :param zmax: height of both domains (m), scalar.
    :param alpha: angle of the upward-sloping part (degree), scalar.
    """
    dx = 0.25
    freq = 340.
    h_s = 2.
    l1max = 60.
    l2max = 40.
    zmax = 40.
    alpha = 20.
    main_analytic_slope(freq, h_s, l1max, l2max, zmax, alpha, dx)
