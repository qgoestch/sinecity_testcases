# -*- coding: utf-8 -*-
##
# \file     plot_errors_norms_cfa2018.py
# \title    Errors and norms for each case: _cfa2018 only 2-norm and max-norm
# \author   Pierre Chobeau
# \version  0.1
# \date     2018, 15 Mar.
##
import numpy as np
import matplotlib.ticker
from matplotlib import pyplot as plt
import os
base_path = reduce (lambda l,r: l + os.path.sep + r,
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def plot_error_basic_cfa2018(h_set, one_norm, two_norm, max_norm,
                     ord_acc_one, ord_acc_two, ord_acc_max,
                     case, save_fig):
    """
    Main plot made of 3 subplots that show (1) the avaraged error,
    (2) the two-norm of the error and (3) the max-norm of the error.

    cfa2018: show two- and max-norms only.

    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param avg_error: error averaged over all receivers for
    each spatial step.
    :type avg_error_tlm: 1d-array
    :param two_norm: relative error in the 2-norm for
    each spatial step.
    :type two_norm: 1d-array
    :param max_norm: relative error in the MAX-norm for
    each spatial step.
    :type max_norm: 1d-array
    :param ord_acc: order of accuracy between two consecutive grids in
    the 2-norm.
    :param case: integer that sorts of the saved folders in the results dir.
    :type case: int
    :param save_fig: save or not the figure.
    :type save_fig: bool
    :type ord_acc: 1d-array

    :return: two graphs: the errors and norms,and the order of accuracy.
    """
    print 'Plotting the errors'
    h_th = np.linspace(h_set[0] - 0.001, h_set[-1] + 0.001, 100)
    j = 2

    # =========================================================================
    #   All grids figure
    # =========================================================================
    fig = plt.figure('Errors', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.loglog(h_set, two_norm[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    error_margin = 0.02 * (h_set[j] / h_set[j]) ** 2 * one_norm[j]
    scnd_ord_th = (h_set / h_set[j]) ** 2 * two_norm[j] + error_margin
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * two_norm[j], 'k-', lw=1.5)
    plt.legend(('FD', '2nd order'), fontsize=14)

    # =========================================================================
    #   Linear regression on log log
    # =========================================================================
    coefs = np.polyfit(h_set, two_norm, 1)
    poly = np.poly1d(coefs)
    ys = poly(h_set)

    m, c = np.polyfit(h_set, np.log10(two_norm), 1)  # fit log(y) = m*log(x) + c
    y_fit = np.power(10, m * h_set + c)  # calculate the fitted values of y
    # print m, c
    # plt.plot(h_set, y_fit, 'y--', lw=3)

    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{2}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.001)
    plt.ylim(10 ** -5, 10 ** -1)
    plt.tight_layout()

    ax = fig.add_subplot(122)
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * max_norm[j], 'k-', lw=1.5)
    ax.loglog(h_set, max_norm[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{max}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.001)
    plt.ylim(10 ** -5, 10 ** -1)
    plt.tight_layout()

    if save_fig:
        print 'Saving figures'
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                'results', 'case%i' % case, 'figures')
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        plt.savefig(os.path.join(res_path, 'errors_fd_cfa2018.eps'),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'errors_fd_cfa2018.png'),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'errors_fd_cfa2018.pdf'),
                    transparent=True, bbox_inches='tight', pad_inches=0)

    # =========================================================================
    #   Order of accuracy btw. 2 consecutive points
    # =========================================================================
    fig = plt.figure('Order of accuracy', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.semilogx(h_set[:-1], ord_acc_two[:], 'bd',
                markersize=4, markeredgewidth=1.2, markeredgecolor='b',
                markerfacecolor='None')
    # plt.legend(('TLM', 'FDTD'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)

    ax = fig.add_subplot(122)
    ax.semilogx(h_set[:-1], ord_acc_max[:], 'bd',
                markersize=4, markeredgewidth=1.2, markeredgecolor='b',
                markerfacecolor='None')
    # plt.legend(('TLM', 'FDTD'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)
    plt.tight_layout()

    if save_fig:
        res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                'results', 'case%i' % case, 'figures')
        if not os.path.exists(res_path):
            os.makedirs(res_path)
        plt.savefig(os.path.join(res_path, 'ord_acc_fd_cfa2018.eps'),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'ord_acc_fd_cfa2018.png'),
                    transparent=True, bbox_inches='tight', pad_inches=0)
        plt.savefig(os.path.join(res_path, 'ord_acc_fd_cfa2018.pdf'),
                    transparent=True, bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_errors_norms_cfa2018(h_set, avg_error_tlm, avg_error_fdtd,
                              two_norm_tlm, two_norm_fdtd,
                              max_norm_tlm, max_norm_fdtd,
                              ord_acc_tlm_one, ord_acc_fdtd_one,
                      ord_acc_tlm_two, ord_acc_fdtd_two,
                      ord_acc_tlm_max, ord_acc_fdtd_max, case):
    """
    Main plot made of 3 subplots that show (1) the avaraged error,
    (2) the two-norm of the error and (3) the max-norm of the error.

    cfa2018: show two- and max-norms only.

    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param avg_error_tlm: error averaged over all receivers for the TLM for
    each spatial step.
    :type avg_error_tlm: 1d-array
    :param avg_error_fdtd: error averaged over all receivers for the FDTD for
    each spatial step.
    :type avg_error_fdtd: 1d-array
    :param two_norm_tlm: relative error in the 2-norm for the TLM for
    each spatial step.
    :type two_norm_tlm: 1d-array
    :param two_norm_fdtd: relative error in the 2-norm for the FDTD for
    each spatial step.
    :type two_norm_fdtd: 1d-array
    :param max_norm_tlm: relative error in the MAX-norm for the TLM for
    each spatial step.
    :type max_norm_tlm: 1d-array
    :param max_norm_fdtd: relative error in the MAX-norm for the FDTD for
    each spatial step.
    :type max_norm_fdtd: 1d-array
    :param ord_acc_tlm_two: order of accuracy between two consecutive grids in
    the 2-norm for the TLM.
    :type ord_acc_tlm_two: 1d-array
    :param ord_acc_fdtd_two: order of accuracy between two consecutive grids in
    the 2-norm for the FDTD.
    :type ord_acc_fdtd_two: 1d-array
    :param ord_acc_tlm_max: order of accuracy between two consecutive grids in
    the max-norm for the TLM.
    :type ord_acc_tlm_max: 1d-array
    :param ord_acc_fdtd_max: order of accuracy between two consecutive grids in
    the max-norm for the FDTD.
    :type ord_acc_fdtd_max: 1d-array
    :param case: integer that sorts of the saved folders in the results
    directory.
    :type case: int
    :return: two graphs, first the errors and norms, second the order of
    accuracy for each norm.
    """
    print 'Plotting the errors'
    h_th = np.linspace(h_set[0] - 0.001, h_set[-1] + 0.001, 100)
    j = 1
    fig = plt.figure('Errors', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.loglog(h_set, two_norm_tlm[:], 'rs',
              markersize=7, markeredgewidth=1.8, markeredgecolor='r',
              markerfacecolor='None')
    ax.loglog(h_set, two_norm_fdtd[:], 'go',
              markersize=4, markeredgewidth=1.8, markeredgecolor='g',
              markerfacecolor='None')
    ax.loglog(h_th, (h_th / h_set[j]) ** 1 * two_norm_tlm[j], 'm--', lw=1)
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * two_norm_tlm[j], 'b-', lw=1)
    plt.legend(('TLM', 'FDTD', '1st order', '2nd order'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{2}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.01)
    plt.ylim(10 ** -5, 10 ** -1)
    plt.tight_layout()

    ax = fig.add_subplot(122)
    ax.loglog(h_th, (h_th / h_set[j]) ** 1 * max_norm_tlm[j], 'm--', lw=1)
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * max_norm_tlm[j], 'b-', lw=1)
    ax.loglog(h_set, max_norm_tlm[:], 'rs',
              markersize=7, markeredgewidth=1.8, markeredgecolor='r',
              markerfacecolor='None')
    ax.loglog(h_set, max_norm_fdtd[:], 'go',
              markersize=4, markeredgewidth=1.8, markeredgecolor='g',
              markerfacecolor='None')
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{max}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.01)
    plt.ylim(10 ** -5, 10 ** -1)
    plt.tight_layout()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'errors_cfa2018.eps'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'errors_cfa2018.png'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'errors_cfa2018.pdf'), transparent=True,
                bbox_inches='tight', pad_inches=0)

    # =========================================================================
    #   Order of accuracy btw. 2 consecutive points
    # =========================================================================
    fig = plt.figure('Order of accuracy', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.semilogx(h_set[:-1], ord_acc_tlm_two[:], 'rs',
                markersize=7, markeredgewidth=1.8, markeredgecolor='r',
                markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_fdtd_two[:], 'go',
                markersize=4, markeredgewidth=1.8, markeredgecolor='g',
                markerfacecolor='None')
    plt.legend(('TLM', 'FDTD'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)

    ax = fig.add_subplot(122)
    ax.semilogx(h_set[:-1], ord_acc_tlm_max[:], 'rs',
                markersize=7, markeredgewidth=1.8, markeredgecolor='r',
                markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_fdtd_max[:], 'go',
                markersize=4, markeredgewidth=1.8, markeredgecolor='g',
                markerfacecolor='None')
    # plt.legend(('TLM', 'FDTD'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)
    plt.tight_layout()
    plt.savefig(os.path.join(res_path, 'ord_acc_cfa2018.eps'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'ord_acc_cfa2018.png'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'ord_acc_cfa2018.pdf'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_errors_norms_fd_fdtd_tlm_cfa2018(h_set, one_norm_fd,
                                          one_norm_tlm, one_norm_fdtd,
                                  two_norm_fd, two_norm_tlm, two_norm_fdtd,
                                  max_norm_fd, max_norm_tlm, max_norm_fdtd,
                                  ord_acc_fd_one, ord_acc_tlm_one, ord_acc_fdtd_one,
                                  ord_acc_fd_two, ord_acc_tlm_two, ord_acc_fdtd_two,
                                  ord_acc_fd_max, ord_acc_tlm_max, ord_acc_fdtd_max, case):
    """
    Main plot made of 3 subplots that show (1) the avaraged error,
    (2) the two-norm of the error and (3) the max-norm of the error.

    cfa2018: show two- and max-norms only.

    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param avg_error_tlm: error averaged over all receivers for the TLM for
    each spatial step.
    :type avg_error_tlm: 1d-array
    :param avg_error_fdtd: error averaged over all receivers for the FDTD for
    each spatial step.
    :type avg_error_fdtd: 1d-array
    :param two_norm_tlm: relative error in the 2-norm for the TLM for
    each spatial step.
    :type two_norm_tlm: 1d-array
    :param two_norm_fdtd: relative error in the 2-norm for the FDTD for
    each spatial step.
    :type two_norm_fdtd: 1d-array
    :param max_norm_tlm: relative error in the MAX-norm for the TLM for
    each spatial step.
    :type max_norm_tlm: 1d-array
    :param max_norm_fdtd: relative error in the MAX-norm for the FDTD for
    each spatial step.
    :type max_norm_fdtd: 1d-array
    :param ord_acc_tlm_two: order of accuracy between two consecutive grids in
    the 2-norm for the TLM.
    :type ord_acc_tlm_two: 1d-array
    :param ord_acc_fdtd_two: order of accuracy between two consecutive grids in
    the 2-norm for the FDTD.
    :type ord_acc_fdtd_two: 1d-array
    :param ord_acc_tlm_max: order of accuracy between two consecutive grids in
    the max-norm for the TLM.
    :type ord_acc_tlm_max: 1d-array
    :param ord_acc_fdtd_max: order of accuracy between two consecutive grids in
    the max-norm for the FDTD.
    :type ord_acc_fdtd_max: 1d-array
    :param case: integer that sorts of the saved folders in the results
    directory.
    :type case: int
    :return: two graphs, first the errors and norms, second the order of
    accuracy for each norm.
    """
    print 'Plotting the errors'
    h_th = np.linspace(h_set[0] - 0.001, h_set[-1] + 0.001, 100)
    j = 1
    fig = plt.figure('Errors', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.loglog(h_set, two_norm_fd[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    ax.loglog(h_set, two_norm_fdtd[:], 'go',
              markersize=5, markeredgewidth=1.8, markeredgecolor='g',
              markerfacecolor='None')
    ax.loglog(h_set, two_norm_tlm[:], 'rs',
              markersize=5, markeredgewidth=1.8, markeredgecolor='r',
              markerfacecolor='None')
    ax.loglog(h_th, (h_th / h_set[j]) ** 1 * two_norm_tlm[j], 'm--', lw=1)
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * two_norm_tlm[j], 'b-', lw=1)
    plt.legend(('FD', 'FDTD', 'TLM', '1st order', '2nd order'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{2}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.01)
    plt.ylim(10 ** -8, 10 ** -0)
    plt.tight_layout()

    ax = fig.add_subplot(122)
    ax.loglog(h_th, (h_th / h_set[j]) ** 1 * max_norm_tlm[j], 'm--', lw=1)
    ax.loglog(h_th, (h_th / h_set[j]) ** 2 * max_norm_tlm[j], 'b-', lw=1)
    ax.loglog(h_set, max_norm_fd[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    ax.loglog(h_set, max_norm_fdtd[:], 'go',
              markersize=5, markeredgewidth=1.8, markeredgecolor='g',
              markerfacecolor='None')
    ax.loglog(h_set, max_norm_tlm[:], 'rs',
              markersize=5, markeredgewidth=1.8, markeredgecolor='r',
              markerfacecolor='None')
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel(r'$||error||_{max}$', fontsize=12)
    plt.xlim(h_set[0] - 0.001, h_set[-1] + 0.01)
    plt.ylim(10 ** -8, 10 ** -0)
    plt.tight_layout()

    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results', 'case%i' % case, 'figures')
    if not os.path.exists(res_path):
        os.makedirs(res_path)
    plt.savefig(os.path.join(res_path, 'errors_3_cfa2018.eps'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'errors_3_cfa2018.png'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'errors_3_cfa2018.pdf'), transparent=True,
                bbox_inches='tight', pad_inches=0)

    # =========================================================================
    #   Order of accuracy btw. 2 consecutive points
    # =========================================================================
    fig = plt.figure('Order of accuracy', figsize=(9., 4.))
    ax = fig.add_subplot(121)
    ax.semilogx(h_set[:-1], ord_acc_fd_two[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_fdtd_two[:], 'go',
                markersize=5, markeredgewidth=1.8, markeredgecolor='g',
                markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_tlm_two[:], 'rs',
                markersize=5, markeredgewidth=1.8, markeredgecolor='r',
                markerfacecolor='None')
    plt.legend(('FD', 'FDTD', 'TLM'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)

    ax = fig.add_subplot(122)
    ax.semilogx(h_set[:-1], ord_acc_fd_max[:], 'bd',
              markersize=4, markeredgewidth=1.2, markeredgecolor='b',
              markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_fdtd_max[:], 'go',
                markersize=5, markeredgewidth=1.8, markeredgecolor='g',
                markerfacecolor='None')
    ax.semilogx(h_set[:-1], ord_acc_tlm_max[:], 'rs',
                markersize=5, markeredgewidth=1.8, markeredgecolor='r',
                markerfacecolor='None')
    # plt.legend(('TLM', 'FDTD'))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    ax.grid(True, which="both", ls=":")
    plt.xlabel('$h$ (m)', fontsize=12)
    plt.ylabel('Obs. order of accuracy', fontsize=12)
    plt.ylim(0, 4)
    plt.tight_layout()
    plt.savefig(os.path.join(res_path, 'ord_acc_3_cfa2018.eps'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'ord_acc_3_cfa2018.png'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.savefig(os.path.join(res_path, 'ord_acc_3_cfa2018.pdf'), transparent=True,
                bbox_inches='tight', pad_inches=0)
    plt.show()
