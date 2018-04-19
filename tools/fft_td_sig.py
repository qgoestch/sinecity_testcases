# -*- coding: utf-8 -*-
##
# \file     fft_td_sig.py
# \title    Calculation of the FFT for the time domain pressure signals.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 14 Feb.
##
import numpy as np
import os
base_path = reduce (lambda l,r: l + os.path.sep + r, 
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )

def post_proc_fft(num_method,case):
    """
     Basic FFT that can be used for time-domain pressures.
    
    :param  num_method   name of the numerical method, string.
    :param      t   discrete time sequence (s), 1d array.
    :param      Ts  time step (s), scalar.
    :param      p_t total acoustic pressure (Pa), 1d array.
    :param      p_f free-field acoustic pressure (Pa), 1d array.
    :param      Fs  sampling frequency (Hz), scalar.
    :param      n   length of the time sequence, scalar.
    :param      f   frequency sequence (s), 1d array.
    :param      w   hanning window, 1d array.
    :param      P_T total acoustic pressure as function of frequency, 1d array.
    :param      P_F free-field acoustic pressure as function of frequency, 1d array.
    
    :return     Return the pressures as a function of the frequency.
    """
#==============================================================================
#   Load the saved pressures and additional variables
#==============================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results','case%i' %case,num_method )
    t       = np.load(os.path.join(res_path,'t.npy'))
    Ts      = np.load(os.path.join(res_path,'Ts.npy'))
    p_t     = np.load(os.path.join(res_path,'p_t.npy'))
    p_f     = np.load(os.path.join(res_path,'p_f.npy'))

#==============================================================================
#   Set the parameters for the FFT
#==============================================================================
    Fs      = 1./Ts
    t       = t[:-1];       n=len(t)
    f       = Fs/1.0*np.arange(0,n)/n
    w       = np.hanning(len(p_t))
    P_T     = np.fft.fft(p_t*w)
    P_F     = np.fft.fft(p_f*w)
    
    return t, p_t, w, P_T, P_F, f
    
def fft_conv(num_method,case,h_num):
    """
     Basic FFT that can be used for time-domain pressures.
    
    :param  num_method   name of the numerical method, string.
    :param  h_num   number of the grid size, scalar.
    :param      t   discrete time sequence (s), 1d array.
    :param      Ts  time step (s), scalar.
    :param      p_t total acoustic pressure (Pa), 1d array.
    :param      p_f free-field acoustic pressure (Pa), 1d array.
    :param      Fs  sampling frequency (Hz), scalar.
    :param      n   length of the time sequence, scalar.
    :param      f   frequency sequence (s), 1d array.
    :param      w   hanning window, 1d array.
    :param      P_T total acoustic pressure as function of frequency, 1d array.
    :param      P_F free-field acoustic pressure as function of frequency, 1d array.
    
    :return     Return the pressures as a function of the frequency.
    """
#==============================================================================
#   Load the saved pressures and additional variables
#==============================================================================
    res_path = os.path.join(base_path.rsplit(os.sep, 1)[0],
                            'results','case%i' %case,num_method )
    t       = np.load(os.path.join(res_path,'t_%i.npy' %(h_num)))
    Ts      = np.load(os.path.join(res_path,'Ts_%i.npy' %(h_num)))
    p_t     = np.load(os.path.join(res_path,'p_t_%i.npy' %(h_num)))
    p_f     = np.load(os.path.join(res_path,'p_f_%i.npy' %(h_num)))

#==============================================================================
#   Set the parameters for the FFT
#==============================================================================
    Fs      = 1./Ts
    t       = t[:-1];       n=len(t)
    f       = Fs/1.0*np.arange(0,n)/n
    w       = np.hanning(len(p_t))
    P_T     = np.fft.fft(p_t*w)
    P_F     = np.fft.fft(p_f*w)
    
    return P_T, P_F, f