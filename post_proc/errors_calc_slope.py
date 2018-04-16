# -*- coding: utf-8 -*-
##
# \file     errors_calc_slope.py
# \title    Comparison of the spectrum and convergence plot.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 02 Oct.
##
import numpy as np
import os
import site
from plot_polar_diag_scat import polar_plots

base_path = reduce(lambda l, r: l + os.path.sep + r,
                   os.path.dirname(os.path.realpath(__file__)).split(os.path.sep))

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0], 'num_methods', 'analytic')
site.addsitedir(analytic_path)

data_plotting_path = os.path.join(base_path.rsplit(os.sep, 2)[0], 'data_plotting')
site.addsitedir(data_plotting_path)
from plot_pres_slope_norot import plot_pres_slope_without_rotation

from analytic_up_slop import main_analytic_slope
from fft_td_sig_arrays import fft_conv_arrays, basic_fft


def error_calc(h_set, h_s, alpha, case, method_name):
    """
    
    :param h_set: spatial step sequence (m).
    :type h_set: list of floats
    :param h_s: source heights (m).
    :type h_s: float
    :param alpha: angle of the upward-sloping part (rad).
    :type alpha: float
    :param case: integer that sorts of the saved folders in the results directory.
    :type case: int
    :param method_name: name of the numerical method of interest.
    :type method_name: str

    :param  res_path:    path where to fetch the results for each method, string chain.
    :param  p_t_fdtd:    total pressure in presence of the circular obstacle, 3d list [angle,distance,time].
    :param  p_f_fdtd:    free-field pressure for the fdtd, 3d list [angle,distance,time].
    :param  f_fdtd:      frequency sequence for the fdtd, 1d list of floats (Hz).
    :param  t_fdtd:      time sequence for the fdtd, 1d list of floats (s).
    :param  Ts_fdtd:     time step in the fdtd, float (s).
    :param  P_fdtd:      Fourier transform of the scattered pressure for the FDTD, 3d list [angle,distance,frequency].
    :param  p_t_tlm:     total pressure in presence of the circular obstacle, 3d list [angle,distance,time].
    :param  p_f_tlm:     free-field pressure for the tlm, 3d list [angle,distance,time].
    :param  f_tlm:       frequency sequence for the tlm, 1d list.
    :param  t_tlm:       time sequence for the tlm, 1d list of floats (s).
    :param  Ts_tlm:      time step in the tlm, float (s).
    :param  P_tlm:       Fourier transform of the scattered pressure for the TLM, 3d list [angle,distance,frequency].

    :param  freq_idx:        chosen frequency index within the range(max_freq_idx), scalar.

    :return: The log-log plot of the two-norm and max-norm, and the respective orders of accuracy for the chosen frequency.
    """
    al_idx = 1
    if method_name == 'fdtd':
        for h_num, h in enumerate(h_set[0:1]):
            # h_num = h_num + 4
            print 'h = %f ; %i/%i' % (h,h_num,len(h_set)-1)

            #   Load the numerical results ---FDTD--- and calculation of the FFTs
            res_path_fdtd = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                         'results', 'case%i' % case, 'fdtd')
            Ts_fdtd = np.load(os.path.join(res_path_fdtd, 'Ts_%i.npy' % h_num))
            t_fdtd = np.load(os.path.join(res_path_fdtd, 't_%i.npy' % h_num))

            # load the compressed array as a numpy.lib.npyio.NpzFile object
            print 'ready to load the *.npz'
            p_fdtd = np.load(os.path.join(res_path_fdtd, 'p_%i_%ideg.npz' % (h_num, int(alpha[al_idx]))))
            print 'p_fdtd.npz loaded'

            # contains a list of the stored arrays in the format '<name>.npy'
            namelist_fdtd = p_fdtd.zip.namelist()
            print namelist_fdtd, namelist_fdtd[1], namelist_fdtd[0]
            os.chdir('../results/case5/fdtd')

            # ONLY TO BE RUN ONE TIME PER GEOMETRY
            # # extract 'p_l*_fdtd.npy' into the current directory
            if not os.path.exists('p_l1.npy') and os.path.exists('p_l2.npy'):
                print 'ready to extraxt the p_l*_fdtd.npy, SLOW!!!'
                p_fdtd.zip.extract(namelist_fdtd[1])
                print 'p_l1_fdtd extracted'
                p_fdtd.zip.extract(namelist_fdtd[0])
                print 'p_l2_fdtd extracted'

            # now we can open the array as a memmap
            print 'ready to load, fast'
            p_l1_fdtd_memmap = np.load(namelist_fdtd[1], mmap_mode='r+')
            print 'l1 loaded'
            p_l2_fdtd_memmap = np.load(namelist_fdtd[0], mmap_mode='r+')
            print 'l2 loaded'

            n_idx_step = 5 # how much nodes you want to skip in order to fasten the FFTs
            print 'Calculation of the FFTs, SLOW if n_idx_step < 10 !!!'
            P_l1_fdtd, f_fdtd = basic_fft(p_l1_fdtd_memmap[::n_idx_step, ::n_idx_step, :], t_fdtd, Ts_fdtd, 2, 'hanning')
            P_l2_fdtd, f_fdtd = basic_fft(p_l2_fdtd_memmap[::n_idx_step, ::n_idx_step, :], t_fdtd, Ts_fdtd, 2, 'hanning')
            print 'FFTs done'

        f_idx = 136
        print f_fdtd[f_idx]
        plot_pres_slope_without_rotation(P_l1_fdtd[:, :, f_idx], P_l2_fdtd[:, :, f_idx], int(alpha[al_idx]),
                                         h_s, f_fdtd[f_idx], h, n_idx_step, num_meth='fdtd')

    if method_name == 'tlm':
        for h_num, h in enumerate(h_set[0:1]):
            #   Load the numerical results FOR THE ---TLM--- and calculation of the FFTs
            res_path_tlm = os.path.join(base_path.rsplit(os.sep, 1)[0],
                                        'results', 'case%i' % case, 'tlm')
            Ts_tlm = np.load(os.path.join(res_path_tlm, 'Ts_%i.npy' % h_num))
            t_tlm = np.load(os.path.join(res_path_tlm, 't_%i.npy' % h_num))

            # load the compressed array as a numpy.lib.npyio.NpzFile object
            print 'ready to load the *.npz'
            p_tlm = np.load(os.path.join(res_path_tlm, 'p_%i_%ideg.npz' % (h_num, int(alpha[al_idx]))))
            print 'p_tlm.npz loaded'

            # contains a list of the stored arrays in the format '<name>.npy'
            namelist_tlm = p_tlm.zip.namelist()
            print namelist_tlm, namelist_tlm[1], namelist_tlm[0]
            os.chdir('../results/case5/tlm')

            # ONLY TO BE RUN ONE TIME PER GEOMETRY
            # extract 'p_l*_tlm.npy' into the current directory
            if not os.path.exists('p_l1.npy') and os.path.exists('p_l2.npy'):
                print 'ready to extraxt the p_l*_tlm.npy, SLOW!!!'
                p_tlm.zip.extract(namelist_tlm[1])
                print 'p_l1_tlm extracted'
                p_tlm.zip.extract(namelist_tlm[0])
                print 'p_l2_tlm extracted'

            # now we can open the array as a memmap
            print 'ready to load, fast'
            p_l1_tlm_memmap = np.load(namelist_tlm[1], mmap_mode='r+')
            print 'l1 loaded'
            p_l2_tlm_memmap = np.load(namelist_tlm[0], mmap_mode='r+')
            print 'l2 loaded'

            print 'Calculation of the FFTs, SLOW if n_idx_step < 10 !!!'
            n_idx_step = 5  # how much nodes you want to skip in order to fasten the FFTs
            P_l1_tlm, f_tlm = basic_fft(p_l1_tlm_memmap[::n_idx_step, ::n_idx_step, :], t_tlm, Ts_tlm, 2, 'hanning')
            P_l2_tlm, f_tlm = basic_fft(p_l2_tlm_memmap[::n_idx_step, ::n_idx_step, :], t_tlm, Ts_tlm, 2, 'hanning')
            print 'FFTs done'

        f_idx = 136
        plot_pres_slope_without_rotation(P_l1_tlm[:, :, f_idx], P_l2_tlm[:, :, f_idx], int(alpha[al_idx]),
                                         h_s, f_tlm[f_idx], h, n_idx_step, num_meth='tlm')

    # P_l1_fdtd_an, P_l2_fdtd_an = main_analytic_slope(f_fdtd[f_idx], h_s, l1max, l2max,
    #                                                  hmax, alpha[al_idx], h)