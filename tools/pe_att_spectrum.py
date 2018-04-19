# -*- coding: utf-8 -*-
##
# \file     pe_att_spectrum.py
# \title    Parabolic equation attenuation spectrum comparison.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 21 Nov.
##

import numpy as np
import os
import site
from matplotlib import pyplot as plt

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                             'num_methods', 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_ground_0

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from miki_imp_model import Miki


def attenuation_spectrum(rho, c, h_s, h_r, d_sr, freq, sigma, case):
    """
    Show the attenuation spectrums along the frequency range defined in the
    main module for the PE simulations. The PE results are compared to
    the analytic solution that is calculated in the present function
    using the simulations parameters.

    :param rho: air density (kg.m-3).
    :type rho: float
    :param c: sound speed (m.s-1).
    :type c: float
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param freq: frequency of the simulation (Hz).
    :type freq: 1darray of floats
    :param sigma: specific airflow resistivity (kNm-4s==CGS).
    :type sigma: list of int
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :return: plots of the attenuation spectrums
    """
    #   Load the numerical results
    sig_idx = 0
    res_path_pe = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                'results', 'case%i' % case, 'pe')
    p_t_pe = [np.load(os.path.join(res_path_pe, 'p_t_%iHz_%icgs.npy'
                              % (f, int(sigma[sig_idx])))) for f in freq]
    p_f_pe = [np.load(os.path.join(res_path_pe, 'p_f_%iHz.npy'
                              % f)) for f in freq]
    pdb_pe = [10. * np.log10(p_t_pe[f_idx][1]**2 / p_t_pe[f_idx][0]**2)
              for f_idx in range(len(freq))]
    # pdb_pe = [10. * np.log10(p[f_idx][1]**2 / (2. * 10**-5))
    #           for f_idx in range(len(freq))]

    pdb_pe_r0 = [10. * np.log10(p_t_pe[f_idx][0]**2 / p_f_pe[f_idx][0]**2)
              for f_idx in range(len(freq))]
    pdb_pe_r1 = [10. * np.log10(p_t_pe[f_idx][1]**2 / p_f_pe[f_idx][1]**2)
              for f_idx in range(len(freq))]
    pdb_pe_r2 = [10. * np.log10(p_t_pe[f_idx][2]**2 / p_f_pe[f_idx][2]**2)
              for f_idx in range(len(freq))]


    #   Calculation of the analytic solutions
    omega = 2. * np.pi * freq
    k_f = omega / c
    Zg, k_m = Miki(-1, freq, sigma[sig_idx], rho, c)
    p_f_an_0, p_t_an_0 = analytic_solution_ground_0(d_sr[0], h_s, h_r,
                                                    Zg / (rho * c), k_f)
    p_f_an_1, p_t_an_1 = analytic_solution_ground_0(d_sr[1], h_s, h_r,
                                                    Zg / (rho * c), k_f)
    p_f_an_2, p_t_an_2 = analytic_solution_ground_0(d_sr[2], h_s, h_r,
                                                    Zg / (rho * c), k_f)
    pdb_an = [10. * np.log10(p_t_an_1[f_idx] ** 2 / p_t_an_0[f_idx] ** 2)
              for f_idx in range(len(freq))]
    # pdb_an = [10. * np.log10(p_t_an_1[f_idx] ** 2 / (2. * 10**-5))
    #           for f_idx in range(len(freq))]

    pdb_an_r0 = [10. * np.log10(p_t_an_0[f_idx] ** 2 / p_f_an_0[f_idx] ** 2)
              for f_idx in range(len(freq))]
    pdb_an_r1 = [10. * np.log10(p_t_an_1[f_idx] ** 2 / p_f_an_1[f_idx] ** 2)
              for f_idx in range(len(freq))]
    pdb_an_r2 = [10. * np.log10(p_t_an_2[f_idx] ** 2 / p_f_an_2[f_idx] ** 2)
              for f_idx in range(len(freq))]

    fig = plt.figure('Att. spect. rel. to r0')
    ax = fig.add_subplot(111)
    plt.semilogx(freq, pdb_pe, 'k*')
    plt.semilogx(freq, pdb_an, 'y--')
    plt.legend(('PE', 'Analytic'), fontsize=12)

    fig = plt.figure('Att. spect. rel. to free field')
    ax = fig.add_subplot(311)
    plt.semilogx(freq, pdb_pe_r0, 'k*')
    plt.semilogx(freq, pdb_an_r0, 'y--')
    plt.legend(('PE', 'Analytic'), fontsize=12)

    ax = fig.add_subplot(312)
    plt.semilogx(freq, pdb_pe_r1, 'k*')
    plt.semilogx(freq, pdb_an_r1, 'y--')

    ax = fig.add_subplot(313)
    plt.semilogx(freq, pdb_pe_r2, 'k*')
    plt.semilogx(freq, pdb_an_r2, 'y--')
    plt.show()
