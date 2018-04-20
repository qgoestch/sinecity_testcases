# -*- coding: utf-8 -*-
##
# \file     case3_pe_ground_imp.py
# \title    Study of an acoustic impulse reflected by an impedance ground.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 21 Nov.
##
import numpy as np
import os
import site

pe_path = os.path.join(os.getcwd().rsplit(os.sep, 1)[0], 'num_methods', 'pe')
site.addsitedir(pe_path)
from init_pe_ground import pe_init_impgr

tools_path = os.path.join(os.getcwd().rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from pe_att_spectrum import attenuation_spectrum


def main(d_sr, h_s, h_r, sigma, freq):
    """
    General script that launch the parabolic equation simulations. The main
    simulation parameters are defined from here, then send to the
    initialization module *init_pe_ground.py*.
    The results postprocessing (errors, spectrum...) may also be launched
    from here.

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param sigma: specific airflow resistivity (kNm-4s==kPam-2s==CGS).
    :type sigma: list of int
    :param freq: frequency of the simulation (Hz).
    :type freq: 1darray of floats
    """
    rho = 1.2   # density (kg.m-3)
    c = 340.00  # sound speed (m.s-1)
    lamb = c / freq     # wavelength (m) - list of floats
    delta_x = lamb / 5.  # 0.017   # spatial step (m) - list of floats
    # delta_x = round(340. / 545.56 / 8. +
    #                 abs(340/545.56/8 - 340/545.56/10)/2., 3)
    # delta_x = 0.125
    case = 38   # integer that sorts of the saved folders in the results dir.
    # free_field = False
    disp_inst_p = True     # display the pressure inside the domain, boolean.
    # for f_idx, f in enumerate(freq):
    #     for sig in sigma:
    #         for ff in [True, False]:
    #             pe_init_impgr(delta_x[f_idx], d_sr, h_s, h_r, f, rho, c,
    #                           sig, case, ff, disp_inst_p)
    #             # pe_init_impgr(delta_x[f_idx], d_sr, h_s, h_r, f, rho, c,
    #             #               sig, case, ff, disp_inst_p)

    attenuation_spectrum(rho, c, h_s, h_r, d_sr, freq, sigma, case)


if __name__ == '__main__':
    d_sr = [10., 60., 100.]
    h_s = 2.
    h_r = 2.
    freq = np.logspace(np.log10(100.), np.log10(428.13291827), 13)
    # freq = np.logspace(np.log10(100.), np.log10(545.559), 15)
    # freq = np.logspace(np.log10(100.), np.log10(1000.), 20)
    # freq = np.arange(50., 275., 25)
    sigma = [200, 20000]
    main(d_sr, h_s, h_r, sigma, freq)
