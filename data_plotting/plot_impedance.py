# -*- coding: utf-8 -*-
##
# \file     plot_impedance.py
# \title    Show the real and imaginary parts of the surface impedance.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 17 Oct.
##
import numpy as np
from scipy import special as sp
import os
import site
from matplotlib import pyplot as plt
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )

tools_path = os.path.join(base_path.rsplit(os.sep, 1)[0],'tools')
site.addsitedir(tools_path)
from miki_imp_model import Miki
from get_imped_coefts import get_coefts_Miki

def plot_surface_imp(sigma, rho, c, f, f_max_src):
    """

    :param sigma: specific airflow resistivity, float (kNm-4s==CGS).
    :param rho: air density, float (kg.m-3).
    :param c: sound speed, scalar (m.s-1).
    :param f: frequency sequence, 1d list of floats (Hz).
    :param f_max_src: maximal frequency, float (Hz)
    :return: the real and imaginary parts of the surface impedance.
    """

#==============================================================================
#   Check the impedance coefficients (btw. Miki and the coef.)
#==============================================================================
    omega       = 2.*np.pi*f
    k           = omega/c
    Zg, k_miki  = Miki(-1,f,sigma,rho,c)
    K           = 6
    a_k, gamma_k, a_k_ncor= get_coefts_Miki(K, sigma)
    am          = 5.50
    bm          = -0.632
    mu          = (am/((2.*np.pi*sigma)**bm))/np.sin(((bm +1.)*np.pi)/2.)
    sum_k       = np.zeros((len(omega)),dtype=np.complex128)
    for n in range(K):
        sum_k   += a_k_ncor[n] / (gamma_k[n]- 1j*omega)  # scalar

    Zomega      = np.zeros((len(omega)),dtype=np.complex128)
    Zomega      = rho*c*(1. + (mu/sp.gamma(-bm))*sum_k)

    plt.figure('Surface impedance')
    plt.semilogx(f, np.real(Zg/(rho*c)), 'k-', lw=2)
    plt.semilogx(f, np.imag(Zg/(rho*c)), 'g-', lw=2)
    plt.semilogx(f, np.real(Zomega/(rho*c)), 'ko', markersize=7, markeredgewidth=2, markeredgecolor='k',
                 markerfacecolor='None')
    plt.semilogx(f, np.imag(Zomega/(rho*c)), 'gs', markersize=5, markeredgewidth=2, markeredgecolor='g',
                 markerfacecolor='None')
    plt.xlim(40,f_max_src)
    plt.ylim(0,max(np.real(Zomega[1]/(rho*c)),np.imag(Zomega[1]/(rho*c))))
    plt.show()