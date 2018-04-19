# -*- coding: utf-8 -*-
##
# \file     errors_calc3_geospr.py
# \title    Calculation of the orders of accuracy as a function of the frequency
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 23 Mar.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(
                       os.path.sep))

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from fft_td_sig_arrays import basic_fft_inner_loading
from obs_ord_acc_freq import obs_ord_acc_geospr

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_obs_ord_acc_freq import plot_obs_ord_acc_freq_geospr


def error_calc3(h_set, case):
    """
    Calculation of the orders of accuracy as a function of frequency from
    the FFTs of the time signals a each receiver.
    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param T:
    :type T:
    :param c: sound speed (m.s-1).
    :type c: float
    :param freq: frequency sequence (Hz).
    :type freq: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    """
    #   Load the numerical results
    P_axi_fdtd, P_dia_fdtd, f_fdtd = zip(*(basic_fft_inner_loading('fdtd', case, i,
                                                                   1, 'Hanning')
                                         for i in range(len(h_set))))
    P_axi_tlm, P_dia_tlm, f_tlm = zip(*(basic_fft_inner_loading('tlm', case, i,
                                                                1, 'Hanning')
                                      for i in range(len(h_set))))
    max_freq_idx = 140
    f_fdtd = f_fdtd[0][:max_freq_idx]
    f_tlm = f_tlm[0][:max_freq_idx]
    print 'FFT done'
    # =========================================================================
    #   Observed order of accuracy as a function of frequency
    # =========================================================================
    p_obs_fdtd_axi = obs_ord_acc_geospr(P_axi_fdtd, f_fdtd)
    p_obs_tlm_axi = obs_ord_acc_geospr(P_axi_tlm, f_tlm)
    p_obs_fdtd_dia = obs_ord_acc_geospr(P_dia_fdtd, f_fdtd)
    p_obs_tlm_dia = obs_ord_acc_geospr(P_dia_tlm, f_tlm)
    print np.shape(p_obs_fdtd_axi)
    plot_obs_ord_acc_freq_geospr(p_obs_fdtd_dia, p_obs_tlm_dia, f_fdtd, case)
    # plot_obs_ord_acc_freq_geospr(P_dia_fdtd[0][:, :max_freq_idx],
    #                              P_dia_tlm[0][:, :max_freq_idx], f_fdtd, case)
    print 'plot done'
