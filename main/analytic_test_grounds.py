# -*- coding: utf-8 -*-
##
# \file     analytic_test_grounds.py
# \title    Test of the analytic solutions - analytic_solutions.py
#           functions: analytic_solution_ground_0,analytic_solution_ground_1,
#           analytic_solution_mixed_ground.
#           Reproduction of the test case shown in [guillaume_jsv2011] for
#           homogenous ground impedance and mixed ground impedances.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 23 Oct.
##
import numpy as np
import os
import site
from scipy import special as sp

analytic_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0],
                             'num_methods','analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_ground_0,\
                                analytic_solution_ground_1,\
                                analytic_solution_mixed_ground
tools_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'tools')
site.addsitedir(tools_path)
from miki_imp_model import Miki


def main(d_sr, h_s, h_r, sigma, f_max_src, d1, d2, h_s2, h_r2, sigma1, sigma2):
    """
    Test of the analytic solutions based on the test case shown in
    **[guillaume_jsv2011, Fig.8 and 9]**, homogenous ground impedance
    **AND** mixed ground impedances.
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param sigma: specific airflow resistivity (kNm-4s==CGS).
    :type sigma: float
    :param f_max_src: the highest frequency of the source signal (Hz).
    :type f_max_src: float
    :param d1: length of the first ground impedance (m).
    :type d1: float
    :param d2: length of the second ground impedance (m).
    :type d2: float
    :param h_s2: height of the source - mixed ground test (m).
    :type h_s2: float
    :param h_r2: height of the receiver - mixed ground test (m).
    :type h_r2: float
    :param sigma1: specific airflow resistivity - ground 1 (kNm-4s==CGS).
    :type sigma1: float
    :param sigma2: specific airflow resistivity - ground 2 (kNm-4s==CGS).
    :type sigma2: float
    :return: plots of the attenuation spectrum for all analytic methods.
    """
    c           = 340.
    rho         = 1.2
    freq        = np.arange(10, f_max_src, 10)
    omega       = 2.*np.pi*freq
    k_freq      = omega/c
    
    # ==========================================================================
    #   IMPEDANCE GROUND [fig.8, guillaume_jsv2011]
    # ==========================================================================
    P_F_analytic_0  = np.zeros((len(freq)),np.complex128)
    P_T_analytic_0  = np.zeros((len(freq)),np.complex128)
    P_F_analytic_1  = np.zeros((len(freq)),np.complex128)
    P_T_analytic_1  = np.zeros((len(freq)),np.complex128)
    P_F_analytic_3  = np.zeros((len(freq)),np.complex128)
    P_T_analytic_3  = np.zeros((len(freq)),np.complex128)

    Zg, k = Miki(-1, freq, sigma, rho, c)
    print 'Calculation in progress... 1/2'
    for i,kfreq in enumerate(k_freq):
        P_F_analytic_0[i], P_T_analytic_0[i] = \
        analytic_solution_ground_0(d_sr,h_s,h_r,Zg[i]/(rho*c),kfreq)
        P_F_analytic_1[i], P_T_analytic_1[i] = \
        analytic_solution_ground_1(d_sr,h_s,h_r,Zg[i]/(rho*c),kfreq)
        P_F_analytic_3[i], P_T_analytic_3[i] = \
        analytic_solution_mixed_ground(10,10,h_s,h_r,Zg[i]/(rho*c),
                                       Zg[i]/(rho*c),kfreq)

    ATT_analytic_0  = 10. * np.log10(np.abs( P_T_analytic_0/P_F_analytic_0)**2)
    ATT_analytic_1  = 10. * np.log10(np.abs( P_T_analytic_1/P_F_analytic_1)**2)
    ATT_analytic_3  = 10. * np.log10(np.abs( P_T_analytic_3/P_F_analytic_3)**2)

    from matplotlib import pyplot as plt
    fig = plt.figure('Impedance ground [fig.8, guillaume_jsv2011]')
    ax  = fig.add_subplot(111)
    ax.plot(freq,ATT_analytic_0, 'r+')
    ax.plot(freq,ATT_analytic_1, 'yx')
    ax.plot(freq,ATT_analytic_3, 'g--')
    ax.set_xticks(np.arange(0,f_max_src,200))
    plt.legend(('Analytic ground 0', 'Analytic ground 1', 'Analytic Mixed'))
    plt.grid()
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Attenuation (dB)")
    plt.xlim(50,f_max_src)
    plt.ylim(-10,6)

    plt.figure('Impedance')
    plt.semilogx(freq, np.real(Zg/(rho*c)), 'k-', lw=2)
    plt.semilogx(freq, np.imag(Zg/(rho*c)), 'g-', lw=2)
    plt.legend((r'Re[$Zg/(\rho*c)$]', r'Im[$Zg/(\rho*c)$]'))
    plt.grid()
    plt.xlabel("Frequency (Hz)")
    plt.ylabel(r'$Zg/(\rho*c)$')
    plt.xlim(10, f_max_src)

    # ==========================================================================
    #   MIXED GROUND [fig.9, guillaume_jsv2011]
    # ==========================================================================
    P_F_analytic = np.zeros(len(freq), np.complex128)
    P_T_analytic = np.zeros(len(freq), np.complex128)

    Zg1, k1 = Miki(-1, freq, sigma1, rho, c)
    Zg2, k2 = Miki(-1, freq, sigma2, rho, c)
    n = 0
    print 'Calculation in progress... 2/2'
    for i, kfreq in enumerate(k_freq[n:]):
        P_F_analytic[i + n], P_T_analytic[i + n] = \
        analytic_solution_mixed_ground(d1, d2, h_s2, h_r2,
                                       Zg1[i + n] / (rho * c),
                                       Zg2[i + n] / (rho * c), kfreq)

    ATT_analytic = 10. * np.log10(np.abs(P_T_analytic / P_F_analytic) ** 2)
    fig = plt.figure('Mixed ground [fig.9, guillaume_jsv2011]')
    ax = fig.add_subplot(111)
    ax.set_xticks(np.arange(0, f_max_src, 200))
    plt.plot(freq, ATT_analytic)
    plt.grid()
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Attenuation (dB)")
    plt.xlim(100, f_max_src)
    plt.ylim(-12, 6)

    plt.figure('Mixed impedance')
    plt.semilogx(freq, np.real(Zg1 / (rho * c)), 'k-', lw=2)
    plt.semilogx(freq, np.imag(Zg1 / (rho * c)), 'g-', lw=2)
    plt.semilogx(freq, np.real(Zg2 / (rho * c)), 'r-', lw=2)
    plt.semilogx(freq, np.imag(Zg2 / (rho * c)), 'm-', lw=2)
    plt.xlim(100, f_max_src)
    plt.show()

if __name__ == '__main__':
    # impedance ground
    d_sr        = 20.
    h_s         = 1.0
    h_r         = 2.0
    sigma       = 50
    f_max_src   = 2000.

    # mixed impedance ground
    d1 = 10.
    d2 = 10.
    h_s2 = 1.
    h_r2 = 2.
    sigma1 = 10
    sigma2 = 1100
    main(d_sr, h_s, h_r, sigma, f_max_src, d1, d2, h_s2, h_r2, sigma1, sigma2)
