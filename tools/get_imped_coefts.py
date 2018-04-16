# -*- coding: utf-8 -*-
##
# \file get_imped_coefts.py
# \title Get coefficients for the recursive convolution method.
# \author Gwenael Guillaume, Nicolas Fortin, Pierre Aumond
# \version 0.1
# \date 2012, 03 august
#
##

import numpy as np
import sys
from math import sin

##
# \fn get_coefts_Miki(K, sigma)
# \brief    Get the coefficients (a_k and gamma_k, residuals and poles of the fraction expansion) 
#           for the given ground parameters in the Miki model.
# \param    K       Order of the partial fraction expansion ;
# \param    sigma   Air flow resistivity (kN.s.m-4).
# \return   The poles and the residuals of the partial faction expansion.
def get_coefts_Miki(K, sigma):
    # Miki impedance model parameters
    am=5.50
    bm=-0.632
    mu=(am/((2.*np.pi*sigma)**bm))/sin(((bm+1.)*np.pi)/2.)
    coefts_dict={}
    coefts_dict['sigma']={}
    if K==6:
        # sigma=10.
        coefts_dict['sigma'][10.]={}       # From article guillaume_jsv2011
        coefts_dict['sigma'][10.]['ak']=[3.45162, -208.80327, 20.64495, 3.53954, 214.40215, 10.34300]
        coefts_dict['sigma'][10.]['lambdak']=[63.25724, 332.11140, 8899.20513, 4.25871, 332.11140, 1660.85954]

        # sigma=20.
        coefts_dict['sigma'][20.]={}       # Validé
        coefts_dict['sigma'][20.]['ak']=[6.816528693543764,76.78398240834527,-1162.146831839459,2183.4952187255335,-1454.6591236065374,400.5321041178154]
        coefts_dict['sigma'][20.]['lambdak']=[4.914502099499762,805.2379210462045,1433.6384964296628,1707.8510040758767,2178.439719528085,3098.857100194904]
        # sigma=50.
        coefts_dict['sigma'][50.]={}       # Validé
        coefts_dict['sigma'][50.]['ak']=[3.44952,-207.87698,20.64699,3.53862,213.47222,10.34131]
        coefts_dict['sigma'][50.]['lambdak']=[63.19759,330.45139,8893.34839,4.25600,330.48282,1658.83025]
        # sigma=60.
        coefts_dict['sigma'][60.]={}       # Validé
        coefts_dict['sigma'][60.]['ak']=[7.504720759233006,18.953907175125348,-100.87252541196693,180.4926503249656,-502.8363377056259,454.21866704591673]
        coefts_dict['sigma'][60.]['lambdak']=[33.83286331755098,679.0797437337246,1479.7848446228409,2051.488869554207,4129.697535722905,4757.87585183627]
        # sigma=80.
        coefts_dict['sigma'][80.]={}       # Validé
        coefts_dict['sigma'][80.]['ak']=[6.283714176198783,10.16737467383077,-2.998197369642349,15.357110599901,-8.40254785293368,47.89699314518404]
        coefts_dict['sigma'][80.]['lambdak']=[0.205917992518911,424.87488208779433,668.6604101324706,2329.8130733780217,5135.941860624568,13515.284123574109]
        # sigma=100.
        coefts_dict['sigma'][100.]={}       # Validé
        coefts_dict['sigma'][100.]['ak']=[7.812193931729983,14.020341833560483,-22.38979709272057,720.4310188656088,-1717.0051863826832,1057.9841838673522]
        coefts_dict['sigma'][100.]['lambdak']=[24.62393143948225,623.5863657378989,1350.2762154733489,3677.2140119387627,4381.08329856621,4985.3801225645275]
        # sigma=200.
        coefts_dict['sigma'][200.]={}       # Validé
        coefts_dict['sigma'][200.]['ak']=[8.31463,256.78346,-370.67054,315.5539,-800.63247,641.07079]
        coefts_dict['sigma'][200.]['lambdak']=[51.08537,991.36901,1111.59341,1807.51167,2885.93181,3268.14348]
        # sigma=300.
        coefts_dict['sigma'][300.]={}       # Validé
        coefts_dict['sigma'][300.]['ak']=[6.37479280139269,8.71608664097787,-415.5025230883408,468.04961907447245,-443.98833789754684,427.02840039146383]
        coefts_dict['sigma'][300.]['lambdak']=[7.89070148160818,384.3421807880428,1438.2339434242,1497.8108983555048,3275.2712415701576,3621.581695884238]
        # sigma=400.
        coefts_dict['sigma'][400.]={}       # Validé
        coefts_dict['sigma'][400.]['ak']=[7.637364963249111,23.50336111285776,-57.25563432492362,253.39282715395524,-493.94314108124547,320.0719668838092]
        coefts_dict['sigma'][400.]['lambdak']=[40.71553866982621,705.1911616745273,1204.2337609158856,2300.094613625184,3213.424781722309,4075.4897980049313]

        coefts_dict['sigma'][490.]={}       # ??? from G. Guillaume, v2, VALID ---> OK
        coefts_dict['sigma'][490.]['ak']=[3.93217, 2.95952, 10.30912, 3.35388, 6.05669, 19.52309]
        coefts_dict['sigma'][490.]['lambdak']=[218.13004, 51.44250, 2839.38884, 3.73547, 779.27735, 11574.74756]

        coefts_dict['sigma'][610.]={}       # ??? from G. Guillaume, v2, VALID ---> OK
        coefts_dict['sigma'][610.]['ak']=[10.58831, 6.15232, 3.70495, 2.71550, 3.18641, 19.75673]     # [2.65312, 10.98489, 2.88270, 3.85431, 6.50341, 19.69437]
        coefts_dict['sigma'][610.]['lambdak']=[2691.88183, 682.11769, 180.78429, 43.74074, 3.27467, 11462.09253]        # [35.04863, 2821.32952, 2.44240, 161.69886, 673.91099, 11941.88069]

        coefts_dict['sigma'][730.]={}       # ??? from G. Guillaume, v2, VALID ---> OK
        coefts_dict['sigma'][730.]['ak']=[2.85687, 3.27865, 6.26111, 3.89074, 19.53409, 10.61415]     # [3.53954, 10.34300, 20.64494, 5.59888, 266.06868, -262.61707]
        coefts_dict['sigma'][730.]['lambdak']=[47.94090, 3.52203, 756.75054, 203.02321, 11872.03603, 2892.58195]        # [4.25871, 1660.85991, 8899.20591, 332.11144, 63.25726, 63.25726]

        # sigma=1100.
        coefts_dict['sigma'][1100.]={}       # Validé
        coefts_dict['sigma'][1100.]['ak']=[8.249597196308462,0.604999243549263,5.78657770709158,9.458053238810747,39.669737056949245,-0.526958354513822]
        coefts_dict['sigma'][1100.]['lambdak']=[60.93420917818878,347.4867349470541,613.7753556020076,2277.4212205652498,12639.69105199048,21879.28455606793]
        # sigma=20000.
        coefts_dict['sigma'][20000.]={}       # Validé
        coefts_dict['sigma'][20000.]['ak']=[0.674623482618556,-3.31744138006447,1.39262011439389,172.443170350865,-214.450105870969,45.8239206927866]
        coefts_dict['sigma'][20000.]['lambdak']=[2.36839619764408,1375.28427802235,1902.76190000849,2000.26500690055,2096.05767353793,2557.53496534917]

    # Correction of poles
    try:
        ak_cor=[(mu/1.42012)*ak_item for ak_item in coefts_dict['sigma'][sigma]['ak']]
        return ak_cor, coefts_dict['sigma'][sigma]['lambdak'], coefts_dict['sigma'][sigma]['ak']
    except KeyError:
        print "The coefficients for an impedance ground for sigma=%.0f cgs are unknown" %(sigma)
        sys.exit() # exit the program

#        elif sigma==80.:
#            ak=[5.529927595219113,6.340797593525356,16.48849236883532,-7.101115871041681,261.1395635263565,-262.2757486761956]
#            lambdak=[6.549074941276979,274.13926422860027,1572.9764533777593,1871.0901904259986,16805.216611106036,22288.699539198005]

##
# \fn get_coefts_Miki_thickness(K, sigma, e)
# \brief Get the coefficients for the given ground parameters in the Miki model
#        with the thickness correction.
# \param K : The order of the partial fraction expansion ;
# \param sigma : The air flow resistivity (kN.s.m-4) ;
# \param e : The material thickness (m).
# \return The poles and the residues.
def get_coefts_Miki_thickness(K, sigma, e):
    # Miki impedance model parameters
    am=5.50
    bm=-0.632
    mu=(am/((2.*np.pi*sigma)**bm))/sin(((bm+1.)*np.pi)/2.)
    coefts_dict={}
    coefts_dict['sigma_e']={}
    if K==5:
        # sigma=200. and e=0.1m
        coefts_dict['sigma_e'][200.,0.1]={}       # Validé
        coefts_dict['sigma_e'][200.,0.1]['ak']=[13342285.1185,-2786526.1138,147587.3634,1964760.0743,14053738.6558]
        coefts_dict['sigma_e'][200.,0.1]['lambdak']=[10.9373,64.4832,1293.9938,22605.2870,22624.8936]
    # Correction of poles
    try:
        ak_cor=[(mu/1.42012)*ak_item for ak_item in coefts_dict['sigma_e'][sigma,e]['ak']]
        return ak_cor,coefts_dict['sigma_e'][sigma,e]['lambdak']
    except KeyError:
        print "The coefficients for an impedance ground as sigma=%.0f cgs are unknown" %(sigma)
        sys.exit() # exit the program





