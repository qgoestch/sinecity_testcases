# -*- coding: utf-8 -*-
##
# \file     obs_ord_acc_freq.py
# \title    Observed order of accuracy as a function of frequency.
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 21 Mar.
##
import numpy as np

def obs_ord_acc(p, h, freq):
    """
    Calculation of the observed order of accuracy as a function of frequency
    from the pressure obtained with 3 grid sizes.
    :param p: pressure from the FDTD method (Pa).
    :type p: 2d-array of floats
    :param h: spatial step sequence (m).
    :type h: list of floats
    :param freq: frequency sequence from the FFT (Hz).
    :type freq: list of floats
    :return: observed order of accuracy as a function of frequency.
    :rtype: list of floats
    """
    p_obs_fdtd = np.zeros((p[0].shape[0], p[0].shape[1],
                           len(freq)))
    for i in range(len(freq)):
        p_obs_fdtd[:, :, i] = np.log(
        (np.abs(p[2][:, :, i]) - np.abs(p[1][:, :, i])) /
        (np.abs(p[1][:, :, i]) - np.abs(p[0][:, :, i]))) / \
        np.log(h[1]/h[0])

    return p_obs_fdtd


def obs_ord_acc_geospr(p, freq):
    """For case 1: geospr.
    Calculation of the observed order of accuracy as a function of frequency
    from the pressure obtained with 3 grid sizes.
    :param p: pressure from the FDTD method (Pa).
    :type p: 2d-array of floats
    :param freq: frequency sequence from the FFT (Hz).
    :type freq: list of floats
    :return: observed order of accuracy as a function of frequency.
    :rtype: list of floats
    """
    p_obs = np.zeros((p[0].shape[0], len(freq)))
    for i in range(len(freq)):
        p_obs[:, i] = \
        (np.abs(p[2][:, i]) - np.abs(p[1][:, i])) / \
        (np.abs(p[1][:, i]) - np.abs(p[0][:, i]))

    return p_obs
