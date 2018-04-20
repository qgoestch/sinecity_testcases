# -*- coding: utf-8 -*-
##
# \file     source_calib.py
# \title    Calibration (normalization) of the Gaussian source signal
#           magnitude in both TLM and FDTD methods.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 29 Jun.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

fdtd_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'fdtd')
site.addsitedir(fdtd_path)
from init_fdtd import fdtd_srl_init

tlm_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'tlm')
site.addsitedir(tlm_path)
from init_tlm import tlm_srl_init

post_proc_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'post_proc')
site.addsitedir(post_proc_path)
from plot_src import plot_src_norm


def main(d_sr, h_s, h_r, sigma, f_max_src):
    """
         Pass the five main input parameters into each numerical method.
      Each method (FDTD or TLM) is launched above a ground, then in free field. The analytic
            solutions are calculated at the last step, in plot_spectrum.py, in which the FFTs are compared.

    :param  d_sr    horizontal distance between the source and the receiver (m), scalar.
    :param  h_s     height of the source (m), scalar.
    :param  h_r     height of the receiver (m), scalar.
    :param  sigma   specific airflow resistivity (kNm-4s==CGS), scalar.
    :param  f_max_src   approximated maximale frequency of the source signal (Hz), scalar.

    :param      c       sound speed, scalar (m.s-1).
    :param      rho     air density, scalar (kg.m-3).
    :param      ref_path_R2   reflected-path length, scalar (m).
    :param      T       simulation duration, scalar (s).
    :param      h_set   spatial step (m).

    :param      case    integer that sorts of the saved folders in the results
                        directory.

    :param      free_field      the domain is enlarged, boolean (True or False).
    :param      disp_inst_p     display the instantaneous pressure, boolean.
    """
    c = 340.
    rho = 1.2
    ref_path_R2 = np.sqrt((d_sr / 2.) ** 2 + h_s ** 2) \
                  + np.sqrt((d_sr / 2.) ** 2 + h_r ** 2)
    T = 2 * (ref_path_R2 / c)
    case = 0
    h_set = 0.001

    fdtd_srl_init(d_sr,h_s,h_r,T,f_max_src,rho,c,sigma,case,True,False)
    tlm_srl_init(d_sr,h_s,h_r,T,f_max_src,rho,c,sigma,case,True,False)

    plot_src_norm(case)

if __name__ == '__main__':
    h_s = 2.0
    d_sr = .01353
    h_r = 2.0
    sigma = 20000.
    f_max_src = 900.

    main(d_sr, h_s, h_r, sigma, f_max_src)