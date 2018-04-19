# -*- coding: utf-8 -*-
##
# \file     plot_error_space.py
# \title    Comparison of the spectrum and convergence plot.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 16 Jun.
##
import numpy as np
from scipy import special as sp
import os
import site
base_path = reduce (lambda l,r: l + os.path.sep + r, 
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )

analytic_path = os.path.join(base_path.rsplit(os.sep, 1)[0],'num_methods','analytic')
site.addsitedir(analytic_path)
from analytic_solutions import analytic_solution_ground_arrays
from fft_td_sig_arrays import fft_conv_arrays
from miki_imp_model import Miki
from error_norm_freq import error, norm, error_n_norm, error_NORMALIZED,\
                            norm_NORMALIZED, pressure_NORMALIZED
def error_space_f(h_set,sigma,rho,c,d_sr,h_s,h_r,f_max_src,case):
    """
      Plot the error and the difference for each method (TLM and FDTD).

    :param  h_step  spatial step (m).
    :param  sigma   specific airflow resistivity (kNm-4s==CGS), scalar.
    :param  rho     air density, scalar (kg.m-3).
    :param  c       sound speed, scalar (m.s-1).
    :param  d_sr    horizontal distance between the source and the receiver (m), scalar.
    :param  h_s     height of the source (m), scalar.
    :param  h_r     height of the receiver (m), scalar.
    :param  f_max_src   approximated maximale frequency of the source signal (Hz), scalar.
    :param  case    integer that sorts of the saved folders in the results
                        directory.

    :param  P_T_fdtd    total pressure in presence of BC for the fdtd, 4d list [h_set][d_sr,h_r,n].
    :param  P_F_fdtd    free-field pressure for the fdtd, 4d list [h_set][d_sr,h_r,n].
    :param  f_fdtd      frequency sequence for the fdtd, 1d list.
    :param  P_T_tlm     total pressure in presence of BC for the tlm, 4d list [h_set][d_sr,h_r,n].
    :param  P_F_tlm     free-field pressure for the tlm, 4d list [h_set][d_sr,h_r,n].
    :param  f_tlm       frequency sequence for the tlm, 1d list.
    
    :param  omega       angular frequency, 1d-list.
    :param  k_f         wave number, 1d-list.
    :param  Zg          surface impedance from Miki's model, 1d-list.
    :param  k_m         wave number from Miki's model, 1d list.
    :param  P_T_an      total pressure in presence of BC - analytic solution, 4d list [h_set][d_sr,h_r,n].
    :param  P_F_an      free-field pressure - analytic solution, 4d list [h_set][d_sr,h_r,n].

    :param  f_end       maximal frequency of study, scalar (int).
    :param  flist       list of the limited frequency sequence, 1d list (int).

    :param  two_norm_tlm    relative error in the 2-norm for the TLM as function
                            of frequency for each spatial step, 1d-array.
    :param  two_norm_fdtd   relative error in the 2-norm for the FDTD as function
                            of frequency for each spatial step, 1d-array.
    :param  two_norm_tlm_MOY    averaged relative error in the 2-norm for the TLM
                                as function of frequency for each spatial step, 1d-array.
    :param  two_norm_fdtd_MOY   averaged relative error in the 2-norm for the FDTD
                                as function of frequency for each spatial step, 1d-array.

    :param  rel_err_tlm     relative error for the TLM as function of frequency
                            for each spatial step, 1d-array.
    :param  rel_err_fdtd    relative error for the FDTD as function of frequency
                            for each spatial step, 1d-array.
    :param  rel_err_tlm_MOY     averaged relative error for the TLM for the coarsest
                                spatial step, 1d-array.
    :param  rel_err_fdtd_MOY    averaged relative error for the FDTD for the coarsest
                                spatial step, 1d-array.
    """
    
#   Load the numerical results
    P_T_fdtd,P_F_fdtd,f_fdtd=zip(*(fft_conv_arrays('fdtd',case,i) for i in range(len(h_set))))
    P_T_tlm ,P_F_tlm ,f_tlm =zip(*(fft_conv_arrays('tlm' ,case,i) for i in range(len(h_set))))
    
#   Calculation of the analytic solutions
    f_end   = 500
    omega   = [2.*np.pi*f_tlm[i] for i in range(len(h_set))]
    k_f     = [omega[i]/c for i in range(len(h_set))]
    Zg, k_m = zip(*(Miki(-1,f_tlm[i],sigma,rho,c) for i in range(len(h_set))))
    print "analytic start"
    P_F_an, P_T_an =zip(*(analytic_solution_ground_arrays(d_sr,h_s,h_r,Zg[i]/(rho*c),k_f[i])
                    for i in range(len(h_set))))
    print "analytic done"
    
#   Errors and norms
    two_norm_tlm_MOY = np.zeros((len(h_set)))
    two_norm_fdtd_MOY= np.zeros((len(h_set)))
    for l in range(len(h_set)):
        flist   = [ int(f_tlm[l][i]) for i in range(1,np.array(P_F_tlm[l]).shape[2]) 
                    if f_tlm[l][i]<f_end]
        rel_err_tlm       = np.zeros((np.array(P_F_an[l]).shape[0],
                                      np.array(P_F_an[l]).shape[1],
                                      len(flist)))
        rel_err_fdtd      = np.zeros((np.array(P_F_an[l]).shape[0],
                                      np.array(P_F_an[l]).shape[1],
                                      len(flist)))
        two_norm_tlm      = np.zeros((len(flist)))
        max_norm_tlm      = np.zeros((len(flist)))
        two_norm_fdtd     = np.zeros((len(flist)))
        max_norm_fdtd     = np.zeros((len(flist)))
        for k in range(len(flist)):
                rel_err_tlm[:,:,k]=\
                error_NORMALIZED(P_T_tlm[l][:,:,k]/P_F_tlm[l][:,:,k],
                                 P_T_an[l][:,:,k]/P_F_an[l][:,:,k])
                two_norm_tlm[k], max_norm_tlm[k] =\
                norm_NORMALIZED(rel_err_tlm[:,:,k],h_set[l])
                rel_err_fdtd[:,:,k]=\
                error_NORMALIZED(P_T_fdtd[l][:,:,k]/P_F_fdtd[l][:,:,k],
                                 P_T_an[l][:,:,k]/P_F_an[l][:,:,k])
                two_norm_fdtd[k], max_norm_fdtd[k] =\
                norm_NORMALIZED(rel_err_fdtd[:,:,k],h_set[l])
        two_norm_tlm    = np.delete(two_norm_tlm, [0, 1])
        two_norm_fdtd   = np.delete(two_norm_fdtd, [0, 1])
        two_norm_tlm_MOY[l] = np.mean(two_norm_tlm)
        two_norm_fdtd_MOY[l]= np.mean(two_norm_fdtd)

    rel_err_tlm_MOY     = np.mean(rel_err_tlm, axis=0)
    rel_err_tlm_MOY     = np.mean(rel_err_tlm_MOY, axis=0)
    rel_err_fdtd_MOY    = np.mean(rel_err_fdtd, axis=0)
    rel_err_fdtd_MOY    = np.mean(rel_err_fdtd_MOY, axis=0)
    
    z_min=0
    z_max=5.

    oo = 1
    pp = int((len(flist)-1)/2.)
    qq = len(flist)-1

    levels = np.arange(z_min,z_max,.01).tolist()
    d_sr_m = [i * .2 for i in range(np.array(P_F_tlm[l]).shape[0])]
    h_r_m  = [i * .2 for i in range(np.array(P_F_tlm[l]).shape[1])]

    X, Y = np.meshgrid(d_sr_m,h_r_m)
    import matplotlib.pyplot as plt
    fig = plt.figure(1)
    ax = fig.add_subplot(311)
    plt.grid()
    plt.contourf(X,Y, np.transpose(rel_err_tlm[:,:,oo]),levels,\
                N=80, norm=None,cmap='magma_r',\
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.1f").set_label(label='rel. error',fontsize=12)
    CS=plt.contour(X,Y, np.transpose(rel_err_tlm[:,:,oo]), levels=np.arange(z_min,z_max,.5), colors='k')
    plt.clabel(CS, fontsize=12, inline=3,fmt = '%.1f',log=True)
    plt.title(r'f = %i Hz' %(flist[oo]), fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    fig.tight_layout()

    ax = fig.add_subplot(312)
    plt.grid()
    plt.contourf(X,Y, np.transpose(rel_err_tlm[:,:,pp]),levels,\
                N=80, norm=None,cmap='magma_r',\
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.1f").set_label(label='rel. error',fontsize=12)
    CS=plt.contour(X,Y, np.transpose(rel_err_tlm[:,:,pp]), levels=np.arange(z_min,z_max,.5), colors='k')
    plt.clabel(CS, fontsize=12, inline=3,fmt = '%.1f',log=True)
    plt.title(r'f = %i Hz' %(flist[pp]), fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    fig.tight_layout()

    ax = fig.add_subplot(313)
    plt.grid()
    plt.contourf(X,Y, np.transpose(rel_err_tlm[:,:,qq]),levels,\
                N=80, norm=None,cmap='magma_r',\
                interpolation='nearest', origin='lower')
    plt.colorbar(format="%.1f").set_label(label='rel. error',fontsize=12)
    CS=plt.contour(X,Y, np.transpose(rel_err_tlm[:,:,qq]), levels=np.arange(z_min,z_max,.5), colors='k')
    plt.clabel(CS, fontsize=12, inline=3,fmt = '%.1f',log=True)
    plt.title(r'f = %i Hz' %(flist[qq]), fontsize=12)
    plt.xlabel('x (m)',fontsize=12)
    plt.ylabel('y (m)',fontsize=12)
    fig.tight_layout()

    plt.show()