# -*- coding: utf-8 -*-
##
# \file     miki_imp_model.py
# \title    Definition of Miki's impedance model that return the characteristic
#           surface impedance used in the analytic models.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 22 Feb.
##
import numpy as np


def Miki(convention, f, sigma_gr, rho_air, c_air):
    """
     Miki's impedance model.
    :param  convention  time convention, '1' for exp(jwt) or '-1' for exp(-jwt)
    :param  f           frequency (Hz), 1d-array.
    :param  e           thickness of the layer (m), scalar.
    :param  sigma_gr    airflow resistivity, !!! in kNm-4s (==CGS) !!!, scalar.
    :param  rho_air     density of air (kg.m-3), scalar.
    :param  c_air       sound speed (m.s-1), scalar.

    :param  K           module d'incompressibilité dynamique du matériau (Pa), scalar.
    :param  sigma_gr    airflow resistivity, !!! in SI: Nm-4s !!!, scalar.
    :param  Zc          surface impedance for semi-inifinite layer (Pasm-1 or kg.m-2.s-1), scalar.
    :param  k           wave number (m-1), scalar.
    :param  Z           impedance for hard-backed layer (thickness correction)  Pasm-1, scalar.
    :param  R           reflection coefficient, scalar.
    :param  alpha       absorbtion coefficient, scalar.
    
    :return Return the layer (ground) surface impedance (without thickness correction).
    """
    Z0      = rho_air*c_air
    sigma_gr_SI   =  sigma_gr * 1000.     # convert the kNm-4s (CGS) into Nm-4s.
    Zc      = Z0*(1.    + 5.5*  (1000.*f/sigma_gr_SI)**-0.632 
                        - 8.43j*(1000.*f/sigma_gr_SI)**-0.632)


    omega   = 2*np.pi*f
    k       = (omega/c_air)*(1. + 7.81*  (1000*f/sigma_gr_SI)**-0.618 
                                - 11.41j*(1000*f/sigma_gr_SI)**-0.618)

##   Densité dynamique
#    rho     = Zc*k/omega
    
##   Module d'incompressibilité
#    K       = omega*Zc/k
    
##   Coefficient d'absorption épaisseur normale
#    Z       = -1j*Zc*np.arctan(k*e)
#    R       = (Z-Z0)/(Z+Z0)
#    alpha   = 1-(abs(R))**2

    if convention == -1:
#        rho = np.conjugate(rho);
#        K   = np.conjugate(K);
        Zc  = np.conjugate(Zc);
#        k   = np.conjugate(k);
#        Z   = np.conjugate(Z);
#        R   = np.conjugate(R);
        
#    np.asarray(Zc)
#    print np.shape(Zc)
    return Zc,k