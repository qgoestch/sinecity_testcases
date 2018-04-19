# -*- coding: utf-8 -*-
##
# \file     source_signals.py
# \title    Definition of the source signals that are used in the initialization
#           of both TLM and FDTD methods.
# \author   Pierre Chobeau
# \version  0.1
# \date     2015, 01 Jan.
##
import numpy as np
import scipy.special as sp


def sine(t, n, freq, pulse_delay):
    """
    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the signal magnitude of a sinusoid at each time iteration
    :rtype: float
    """
    return 1.*np.sin(2*np.pi*freq*t[n])
    

def gauss_1(t, n, freq, pulse_delay):
    """

    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the signal magnitude of a Gaussian pulse at each time iteration
    :rtype: float
    """
    return np.exp(-(np.pi**2) * ((freq/2.)*t[n-pulse_delay]-1)**2)


def gauss_2(t, n, freq, pulse_delay):
    """
    shaeffer, jasa 2014, physically constraint source, Eq.(39)

    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the signal magnitude of a Gaussian pulse at each time iteration
    :rtype: float
    """
    sig=1./3*10**-2
    s=np.sqrt(np.pi/2)*sig*sp.erf(t[n-pulse_delay]/(np.sqrt(2)*sig))
    return s


def ricker(t, n, freq, pulse_delay):
    """

    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the signal magnitude of a Ricker wavelet (also calledmexican hat or
    raised sinus) at each time iteration
    :rtype: float
    """
    pulse_delay=1./freq
    return (1.-(2.*(np.pi**2)*(freq**2)*((t[n]-pulse_delay)**2))) * \
            np.exp(-(np.pi**2)*(freq**2)*((t[n]-pulse_delay)**2))


def dirac(t, n, freq, pulse_delay):
    """

    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the signal magnitude of a Dirac at each time iteration
    :rtype: float
    """
    if n==pulse_delay:
        s=1
    else:
        s=0
    return s


def src_select(src_typ, t, n, freq, pulse_delay):
    """
    Select the source signal using its name and send the parameters.

    :param src_typ: type of source signal
    :type src_typ: string
    :param t: time sequence (s).
    :type t: list of floats
    :param n: time iteration index
    :type n: int
    :param freq: frequency of the sinusoid (Hz)
    :type freq: float
    :param pulse_delay: number of iteration for the delay of the signal defined
    in the init_*.py
    :type pulse_delay: int
    :return: the selected signal magnitude at each time iteration t[n].
    :rtype: float
    """
    if src_typ == 'sine':
        s = sine(t, n, freq, pulse_delay)
    elif src_typ == 'gauss_1':
        s = gauss_1(t, n, freq, pulse_delay)
    elif src_typ == 'gauss_2':
        s = gauss_2(t, n, freq, pulse_delay)
    elif src_typ == 'ricker':
        s = ricker(t, n, freq, pulse_delay)
    elif src_typ == 'dirac':
        s = dirac(t, n, freq, pulse_delay)
    return s