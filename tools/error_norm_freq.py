# -*- coding: utf-8 -*-
##
# \file     error_norm_freq.py
# \title    Calculation of the relative error and the norms.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 13 Apr.
##
import numpy as np
import os
base_path = reduce (lambda l,r: l + os.path.sep + r, 
                    os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep ) )


def error(p_num, p_exa):
    """
     Calculation of the absolute error from the numerical and analytical pressures.

    :param p_num: numerical results (Pa)
    :type p_num: list
    :param p_exa: analytic solution (Pa)
    :type p_exa: list
    :return: The absolute error at a given frequency.
    :rtype: float
    """
    return np.abs(p_num - p_exa)


def one_norm(err, h):
    """
     Calculation of the one-norm from the absolute error.

    :param err: abolute error
    :type err: float
    :param h: spatial step (m)
    :type h: float
    :return: The one-norm of the absolute error at a given frequency.
    :rtype: float
    """
    return np.sum(np.abs(err), dtype=np.float64) * h ** 2


def two_norm(err, h):
    """
     Calculation of the two-norm from the absolute error.

    :param err: abolute error
    :type err: float
    :param h: spatial step (m)
    :type h: float
    :return: The two-norm of the absolute error at a given frequency.
    :rtype: float
    """
    return np.sqrt(np.sum(np.abs(err ** 2), dtype=np.float64) * h ** 2)


def max_norm(err, h):
    """
     Calculation of the max-norm from the absolute error.

    :param err: abolute error
    :type err: float
    :param h: spatial step (m)
    :type h: float
    :return: The max-norm of the absolute error at a given frequency.
    :rtype: float
    """
    return np.sqrt(np.max(err ** 2) * h ** 2)


def error_NORMALIZED(p_num,p_exa):
    """
     Calculation of the relative error as a function of frequency.
    
    :param  p_num   numerical results.
    :param  p_exa   analytic solution.
    :param  f       frequency sequence.
    :param      rel_err relative error as a function of frequency.
    :return     The relative error as a function of frequency.
    """
    # error_L2[k, j, i] = np.abs(np.abs(pfv[k, j, i]) - Mag_Norm * np.abs(pan[k, j, i - 1])) ** 2
    rel_err = np.abs((np.abs(p_num) - np.abs(p_exa))/ np.abs(p_exa))**2
    return rel_err


def norm_NORMALIZED(rel_err,h):
    """
     Calcution of the max-norm and 2-norm as a function of frequency.
    
    :param  rel_err relative error as a function of frequency.
    :param  h       spatial step.
    :param      two_norm  the 2-norm.
    :param      max_norm  the max-norm.
    :return     The 2-norm and the max-norm.
    """
    # L2_Norm[i] = np.sqrt(np.sum(error_L2[:, :, i]) * h ** 2)
    two_norm = np.sqrt(     np.sum(rel_err) *h**2 )
    max_norm = np.sqrt(     np.max(rel_err) *h**2 )
    return two_norm, max_norm


def norm(p_num,p_exa,f,f_idx,rel_err,h,Ts):
    """
     Calcution of the max-norm and 2-norm as a function of frequency.
    
    :param  rel_err relative error as a function of frequency.
    :param  h       spatial step.
    :param  Ts      time step.
    :param      two_norm  the 2-norm.
    :param      max_norm  the max-norm.
    :return     The 2-norm and the max-norm.
    """
    two_norm = np.zeros((rel_err.shape[2]), dtype=np.float64)
    max_norm = np.zeros((rel_err.shape[2]), dtype=np.float64)
    for i in range(rel_err.shape[2]):
        two_norm[i] = np.sqrt(  np.sum(np.abs(rel_err[:,:,i])**2)*h**2)
        max_norm[i] = np.max(   np.abs(rel_err[:,:,i]))*h**2
    return rel_err, two_norm, max_norm


def pressure_NORMALIZED(p_num, norm_ref):
    """
     Calculation of the normalized pressure as a function of frequency.

    :param  p_num   numerical results.
    :param  norm_ref reference pressure for the normalization.
    :param      p_norm  normalized pressure.
    :return     The normalized pressure as a function of frequency.
    """
    p_norm = p_num / norm_ref
    return p_norm

def error_n_norm(p_num,p_exa,f,f_idx,rel_err,h,Ts):
    """
     Calcution of the max-norm and 2-norm.
    
    :param  rel_err relative error as a function of frequency.
    :param  h       spatial step.
    :param  Ts      time step.
    :param      two_norm  the 2-norm.
    :param      max_norm  the max-norm.
    :return     The 2-norm and the max-norm.
    :param  p_num   numerical results.
    :param  p_exa   analytic solution.
    :param  f       frequency sequence.
    :param      rel_err relative error as a function of frequency.
    :return     The relative error as a function of frequency.
    """
    flim_idx = min(p_num.shape[2],p_exa.shape[2]) -1
#==============================================================================
#   Relative error calculation
#==============================================================================
    rel_err_sqr = ( np.abs((p_num[:,:,:flim_idx] - p_exa[:,:,:flim_idx])
                    /p_exa[:,:,:flim_idx]))**2
#    error_L2_NUM[k,j,i] = np.abs(np.abs(pfv[k,j,i]) - Mag_Norm*np.abs(pan[k,j,i-1]))**2
    
    two_norm = np.zeros((flim_idx), dtype=np.float64)
    max_norm = np.zeros((flim_idx), dtype=np.float64)
    for i_idx,i in enumerate(range(flim_idx)):
#        for j in range(rel_err.shape[1]):
#            for k in range(rel_err.shape[0]):
#==============================================================================
#   2-norm: over the receivers (all???)
#==============================================================================
        two_norm[i_idx] = np.sqrt(np.sum(np.abs(rel_err_sqr[:,:,i])) * h**2)
#        L2_Norm[i] = np.sqrt( np.sum(error_L2_NUM[:,:,i]) * h**2 / np.sum(DENOMINATOR[:,:,i]) )

#==============================================================================
#   max-norm
#==============================================================================
        max_norm[i_idx] = np.max(np.abs(rel_err[:,:,i]))*h**2 # frequency only
#    max_norm    = h**2*(1./Ts)*np.max(np.abs(rel_err)**2) # both frequency and point
    
    return rel_err, two_norm, max_norm