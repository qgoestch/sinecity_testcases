# -*- coding: utf-8 -*-
##
# \file     errors_calc2_geospr.py
# \title    Calculation of the errors and norms for the case1:
#           geometrical spreading.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 08 Feb.
##
import numpy as np
from scipy import special as sp
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(
                       os.path.sep))

tools_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'tools')
site.addsitedir(tools_path)
from fft_td_sig_arrays import fft_conv_arrays
from error_norm_freq import error, two_norm, max_norm
import source_signals as src

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_time_signals import plot_ts, plot_ts_basic
from plot_errors_norms import plot_error_basic


def error_calc2(d_sr, h_set, T, c, freq, case):
    """
    :param d_sr: horizontal distances between the source and the receivers,
    list of floats (m).
    :type d_sr: list of floats
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param c: sound speed (m.s-1).
    :type c: float
    :param freq: frequency sequence (Hz).
    :type freq: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    """
    for num_meth in ['fdtd', 'tlm']:
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                                'case%i' % case, num_meth)
        one_norm = np.zeros((len(h_set)))
        two_norm = np.zeros((len(h_set)))
        max_norm = np.zeros((len(h_set)))
        ord_acc_one = np.zeros((len(h_set) - 1))
        ord_acc_two = np.zeros((len(h_set) - 1))
        ord_acc_max = np.zeros((len(h_set) - 1))
        for l in range(len(h_set)):
            p_axial = np.load(os.path.join(res_path, 'p_axial_%i.npy' % l))
            p_diag = np.load(os.path.join(res_path, 'p_diag_%i.npy' % l))

            t = np.load(os.path.join(res_path, 't_%i.npy' % l))
            Ts = np.load(os.path.join(res_path, 'Ts_%i.npy' % l))
            pulse_delay = int(round(1 / 50. / Ts))
            It = range(0, t.shape[0])
            d_idx = 0
            p_an = np.zeros((len(d_sr), It[-1] + 1))
            error = np.zeros((len(d_sr)))
            print d_sr[d_idx]
            for d_idx in range(len(d_sr)):
                for n in It:
                    # 1/14 normalization of the magnitude at the 1st reveiver
                    # using the funtion plot_ts_basic below.
                    p_an[d_idx, n] = 1. / (14. * np.sqrt(d_sr[d_idx])) * \
                                     np.exp(-(np.pi**2) * ((freq/2.) *
                                    (t[n-pulse_delay] - d_sr[d_idx]/c) - 1)**2)
                error[d_idx] = np.mean(p_an[d_idx, :] - p_axial[d_idx, :])
            # plot_ts_basic(t, p_an[d_idx, :], p_diag[d_idx, :])

            one_norm[l] = np.linalg.norm((error) *
                                            h_set[l] ** 2, ord=1)
            two_norm[l] = np.linalg.norm((error) *
                                            h_set[l] ** 2, ord=2)
            max_norm[l] = np.linalg.norm((error) *
                                            h_set[l] ** 2, ord=np.inf)

            for l in range(len(h_set) - 1):
                ord_acc_one[l] = np.log(
                    one_norm[l + 1] / one_norm[l]) / np.log(
                    h_set[l + 1] / h_set[l])
                ord_acc_two[l] = np.log(
                    two_norm[l + 1] / two_norm[l]) / np.log(
                    h_set[l + 1] / h_set[l])
                ord_acc_max[l] = np.log(
                    max_norm[l + 1] / max_norm[l]) / np.log(
                    h_set[l + 1] / h_set[l])

        import os
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                'results', 'case%i'
                                % case, '%s' % num_meth)
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        np.save(os.path.join(res_path, 'one_norm_%s.npy' % num_meth),
                one_norm)
        np.save(os.path.join(res_path, 'two_norm_%s.npy' % num_meth),
                two_norm)
        np.save(os.path.join(res_path, 'max_norm_%s.npy' % num_meth),
                max_norm)
        np.save(os.path.join(res_path, 'ord_acc_one_%s.npy' % num_meth),
                ord_acc_one)
        np.save(os.path.join(res_path, 'ord_acc_two_%s.npy' % num_meth),
                ord_acc_two)
        np.save(os.path.join(res_path, 'ord_acc_max_%s.npy' % num_meth),
                ord_acc_max)
        # plot_error_basic(h_set, one_norm, two_norm, max_norm,
        #                  ord_acc_one, ord_acc_two, ord_acc_max,
        #                  case, True)
