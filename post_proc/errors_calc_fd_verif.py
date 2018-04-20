# -*- coding: utf-8 -*-
##
# \file     errors_calc_fd_verif.py
# \title    Calculation of the errors and norms for the finite difference
#           method applied to the Helmholtz equation.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2018, 02 Feb.
##
import numpy as np
import os
import site

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).
                   split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                             'num_methods', 'analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_ground_arrays

tools_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'tools')
site.addsitedir(tools_path)
from error_norm_freq import error, two_norm, max_norm

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0],
                                  'data_plotting')
site.addsitedir(data_plotting_path)
from plot_errors_norms_cfa2018 import plot_error_basic_cfa2018


def error_calc(h_set, case):
    """
    Calculation of the error for the finite difference applied to the
    Helmholtz equation.

    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int

    :param p_fd: pressure from the finite difference Helmholtz equation (Pa)
    :type p_fd: 2d-array (space) of floats
    :param p_an: pressure from the analytic solution (Pa)
    :type p_an: 2d-array (space) of floats

    :return:
    :rtype:
    """
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'results',
                            'case%i' % case, 'fd')
    # =========================================================================
    #   Calculation of the errors and norms using numpy.linalg.norm
    # Work for the 2-norm, not for the max-norm that is only doable for
    # vectors (not possible for the arrays >= 2D)
    # =========================================================================
    one_norm_fd = np.zeros((len(h_set)))
    two_norm_fd = np.zeros((len(h_set)))
    max_norm_fd = np.zeros((len(h_set)))
    avg_error_fd = np.zeros((len(h_set)))
    ord_acc_one_fd = np.zeros((len(h_set) - 1))
    ord_acc_two_fd = np.zeros((len(h_set) - 1))
    ord_acc_max_fd = np.zeros((len(h_set) - 1))
    for l in range(len(h_set)):
        p_fd = np.load(os.path.join(res_path, 'p_fd_%i.npy' % l))
        p_an = np.load(os.path.join(res_path, 'p_an_%i.npy' % l))
        one_norm_fd[l] = np.linalg.norm((p_fd[1:-1, 1:-1] - p_an[1:-1, 1:-1]) *
                                        h_set[l] ** 2, ord=1)
        two_norm_fd[l] = np.linalg.norm((p_fd[1:-1, 1:-1] - p_an[1:-1, 1:-1]) *
                                        h_set[l] ** 2, ord=2)
        max_norm_fd[l] = np.linalg.norm((p_fd[1:-1, 1:-1] - p_an[1:-1, 1:-1]) *
                                        h_set[l] ** 2, ord=np.inf)

        # e_ij_fd = error(p_fd[1:-1, 1:-1], p_an[1:-1, 1:-1])
        # two_norm_fd[l] = two_norm(e_ij_fd, h_set[l])
        # max_norm_fd[l] = max_norm(e_ij_fd, h_set[l])

    for l in range(len(h_set) - 1):
        ord_acc_one_fd[l] = np.log(
            one_norm_fd[l + 1] / one_norm_fd[l]) / np.log(
            h_set[l + 1] / h_set[l])
        ord_acc_two_fd[l] = np.log(
            two_norm_fd[l + 1] / two_norm_fd[l]) / np.log(
            h_set[l + 1] / h_set[l])
        ord_acc_max_fd[l] = np.log(
            max_norm_fd[l + 1] / max_norm_fd[l]) / np.log(
            h_set[l + 1] / h_set[l])

    num_meth = 'fd'
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i'
                            % case, '%s' % num_meth)
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    np.save(os.path.join(res_path, 'one_norm_%s.npy' % num_meth),
            one_norm_fd)
    np.save(os.path.join(res_path, 'two_norm_%s.npy' % num_meth),
            two_norm_fd)
    np.save(os.path.join(res_path, 'max_norm_%s.npy' % num_meth),
            max_norm_fd)
    np.save(os.path.join(res_path, 'ord_acc_one_%s.npy' % num_meth),
            ord_acc_one_fd)
    np.save(os.path.join(res_path, 'ord_acc_two_%s.npy' % num_meth),
            ord_acc_two_fd)
    np.save(os.path.join(res_path, 'ord_acc_max_%s.npy' % num_meth),
            ord_acc_max_fd)

    plot_error_basic_cfa2018(h_set, one_norm_fd, two_norm_fd, max_norm_fd,
                     ord_acc_one_fd, ord_acc_two_fd, ord_acc_max_fd,
                     case, True)
