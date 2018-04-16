# -*- coding: utf-8 -*-
##
# \file     analytic_solutions.py
# \title    Analytic solutions that gives the pressure above an impedance ground
#           for a point source. The first solution uses erfc function, whereas
#           the second uses the trapezoidal rule to solve the integral.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 24 Oct.
##
import numpy as np
import scipy.special as sp


def analytic_solution_ground_0(d_sr, h_s, h_r, Zg, kfreq):
    """
    Analytic solution for point source emitting plane waves above a ground.
    The ground can be either reflecting or an impedance ground.
    The analytic solution has been published by **[Di & Gilbert in JASA_1993]**.
    See also **[berengier_actaac2003, Eqs.(1-6)]** in which this solution is
    reviewed. Using the error erfc function from scipy.

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param Zg: ground impedance in the frequency domain (Zg(k)).
    :type Zg: 1d numpy array
    :param kfreq: wave number --> omega=2*np.pi*f ; k=omega/c0.
    :type kfreq: 1d numpy array

    :return: Gives both the total pressure and the free field pressure as a
    function of frequency.
    :rtype: 1darray
    """
    R1      = np.sqrt(abs(h_s-h_r)**2 + d_sr**2) # direct path between the source and the receiver, float (m).
    R2      = np.sqrt((h_s+h_r)**2 + d_sr**2) # reflected path assuming a specular reflection, float (m).
    theta   = np.arctan(d_sr/abs(h_s+h_r)) # angle of btween the relfected ptah and the normal, float (rad).
    w       = 0.5*(1.+ 1j)*np.sqrt(kfreq*R2)*(np.cos(theta)+(1./Zg))
    Fw      = 1.+(1j*np.sqrt(np.pi)* w *np.exp(-w**2)* sp.erfc(-1j*w))
    Rp      = (np.cos(theta)-(1/Zg))/(np.cos(theta)+(1./Zg))
    Q       = Rp + (1.-Rp)*Fw
    p_f     = np.exp(1j*kfreq*R1)/(4*np.pi*R1) # free-field pressure, complex (Pa).
    p_t     = np.exp(1j*kfreq*R1)/(4*np.pi*R1)\
            + np.exp(1j*kfreq*R2)/(4*np.pi*R2)*Q # total pressure, complex (Pa).
    return p_f, p_t


def analytic_solution_ground_1(d_sr, h_s, h_r, Zg, kfreq):
    """
     Analytic solution for point source emitting plane waves above a ground.
     The ground can be either reflecting or an impedance ground.
    The analytic solution has been published by **[Di & Gilbert in JASA_1993]**.
    See also [berengier_actaac2003, Eqs.(1-6)] in which this solution is reviewed.
    Using the normalised error erfcx function from scipy.

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param Zg: ground impedance in the frequency domain (Zg(k)).
    :type Zg: 1d numpy array
    :param kfreq: wave number --> omega=2*np.pi*f ; k=omega/c0.
    :type kfreq: 1d numpy array

    :return     Gives both the total pressure and the free field pressure as a function of frequency.
    """
    R1      = np.sqrt((h_r-h_s)**2 + d_sr**2)   # direct path between the source and the receiver, float (m).
    R2      = np.sqrt((h_r+h_s)**2 + d_sr**2)   # reflected path assuming a specular reflection, float (m).
    theta   = np.arctan(d_sr/(h_s+h_r))     # angle of btween the relfected ptah and the normal, float (rad).
    psi     = np.pi/2 - theta
    Rp      = (Zg*np.sin(psi)-1.)/(Zg*np.sin(psi)+1.)
    w       = (1.+1j)/2.*np.sqrt(kfreq*R2)*(np.sin(psi)+(1./Zg))
    Fw      = 1.+1j*np.sqrt(np.pi)*w*sp.erfcx(-1j*w)
    Q       = Rp + (1.-Rp)*Fw
    p_f     = np.exp(1j*kfreq*R1)/(4*np.pi*R1)  # free-field pressure, complex (Pa).
    p_t     = np.exp(1j*kfreq*R1)/(4*np.pi*R1)\
            + np.exp(1j*kfreq*R2)/(4*np.pi*R2)*Q    # total pressure, complex (Pa).
    return p_f,p_t


def analytic_solution_ground_arrays(d_sr, h_s, h_r, Zg, kfreq):
    """
     Copy of analytic_solution_1 but with 1 output.

    :param d_sr: horizontal distances between the source and the receivers (m).
    :type d_sr: list of floats
    :param h_s: height of the source (m).
    :type h_s: float
    :param h_r: height of the receiver (m).
    :type h_r: float
    :param Zg: ground impedance in the frequency domain (Zg(k)).
    :type Zg: 1d numpy array
    :param kfreq: wave number --> omega=2*np.pi*f ; k=omega/c0.
    :type kfreq: 1d numpy array

    :return     Gives both the total pressure and the free field pressure as a function of frequency.
    """
    p_f = np.zeros((len(d_sr), len(h_r), len(kfreq)), np.complex128)
    p_t = np.zeros((len(d_sr), len(h_r), len(kfreq)), np.complex128)
    for k, d in enumerate(d_sr):
        for j, hr in enumerate(h_r):
            for i, kf in enumerate(kfreq):
                R1 = np.sqrt((hr - h_s) ** 2 + d ** 2)  # direct path between the source and the receiver, float (m).
                R2 = np.sqrt((hr + h_s) ** 2 + d ** 2)  # reflected path assuming a specular reflection, float (m).
                theta = np.arctan(d / (h_s + hr))   # angle of btween the relfected ptah and the normal, float (rad).
                psi = np.pi / 2 - theta
                Rp = (Zg[i] * np.sin(psi) - 1.) / (Zg[i] * np.sin(psi) + 1.)
                w = (1. + 1j) / 2. * np.sqrt(kf * R2) * (np.sin(psi) + (1. / Zg[i]))
                #    w       = np.conj(-w*1j)
                Fw = 1. + 1j * np.sqrt(np.pi) * w * sp.erfcx(-1j * w)
                Q       = Rp + (1.-Rp)*Fw
                # Q = 1.
                p_f[k, j, i] = np.exp(1j * kf * R1) / (4 * np.pi * R1)  # free-field pressure, complex (Pa).
                p_t[k, j, i] = np.exp(1j * kf * R1) / (4 * np.pi * R1) \
                               + np.exp(1j * kf * R2) / (4 * np.pi * R2) * Q    # total pressure, complex (Pa).
    return p_f, p_t


def analytic_solution_mixed_ground(d1,d2,h_s,h_r,Zg1,Zg2,kfreq):
    """
     Analytic solution for point source emitting plane waves above a
    **MIXED IMPEDANCE** ground with impedance Z1 and Z2.
     The solution from **[Rasmussen_jsv1982]** is implemented using
    **[berengier_actaac2003, Eqs.(8,9)]** using the trapezoidal rule through the
    numpy library np.trapz.

    :param d1: distance of the first impedance ground - source side (m).
    :type d1: float
    :param d2: distance of the second impedance ground - receiver side (m).
    :type d2: float
    :param h_s: height of the source (m).
    :type d2: float
    :param h_r: height of the receiver (m).
    :type d2: float
    :param Zg1: first ground impedance in the frequency domain as a function
    of the wavenumber - (Pasm-1 or kg.m-2.s-1).
    :type Zg1: 1d numpy array of complexes.
    :param Zg2: second ground impedance in the frequency domain as a function
    of the wavenumber - (Pasm-1 or kg.m-2.s-1).
    :type Zg2: 1d numpy array of complexes.
    :param kfreq: wave number --> omega=2*np.pi*f ; k=omega/c0.
    :type kfreq:  1d numpy array of floats

    :param k1: first ground wave number; 1d numpy array.
    :param k2: second ground wave number; 1d numpy array.

    :param h_z: height of integration at the limit between the 2 grounds; 1d numpy array.
    :param R1: direct path between the source (S) and the limit, scalar (m).
    :param R2: reflected path (ground 1) assuming a specular reflection, scalar (m).
    :param theta1: angle of the specular rectlection (ground 1), scalar (rad).
    :param Rp1: reflection coefficient for the ground 1, scalar.
    :param w1,Fw1: intermediate values for the calculation of the Q values, scalars.
    :param Q1: spherical wave reflection coefficient for the ground 1, scalar.
    :param R3: direct path between the limit and the limit, scalar (m).
    :param R4: reflected path (ground 2) assuming a specular reflection, scalar (m).
    :param theta2: angle of the specular rectlection (ground 2), scalar (rad).
    :param Rp2: reflection coefficient for the ground 2, scalar.
    :param w2,Fw2: intermediate values for the calculation of the Q values, scalars.
    :param Q2: spherical wave reflection coefficient for the ground 2, scalar.
    
    :param p1,p2,p3,p4: pressure for the different paths, floats.
    :param I: integral calculation, float.
    
    :param p_r: reflected pressure (Pa),scalar.
    :param p_t: total pressure field (Pa),scalar.

    :param R_dir: direct path between the source (S) and the receiver (R), scalar (m).
    :param p_f: free-field pressure (Pa),scalar.

    :return the total pressure and the free field pressure as a function of frequency.
    """
    CritF   = 0.07
    lambd   = 2.*np.pi/kfreq
    Zmax    = 7.*(d1+d2)
    dz      = lambd*CritF
    h_z     = np.arange(0,int(Zmax/dz),dz)

    I       = np.zeros(1,np.complex128)
    R1      = np.sqrt((h_z-h_s)**2 + d1**2)
    R2      = np.sqrt((h_z+h_s)**2 + d1**2)
    theta1  = np.arctan(d1/(h_s+h_z))
    psi1    = np.pi/2 + theta1
    Rp1     = (Zg1*np.sin(psi1)-1)/(Zg1*np.sin(psi1)+1)
    w1      = (1.+1j)/2.*np.sqrt(kfreq*R2)*(np.sin(psi1)+(1./Zg1))  
    Fw1     = 1.+1j*np.sqrt(np.pi)*w1*sp.erfcx(-1j*w1)
    Q1      = Rp1 + (1.-Rp1)*Fw1
    
    R3      = np.sqrt((h_z-h_r)**2 + d2**2)
    R4      = np.sqrt(d2**2 + (h_z+h_r)**2 )
    theta2  = np.arctan(d2/(h_r+h_z))
    psi2    = np.pi/2 + theta2
    Rp2     = (Zg2*np.sin(psi2)-1)/(Zg2*np.sin(psi2)+1)
    w2      = (1.+1j)/2.*np.sqrt(kfreq*R4)*(np.sin(psi2)+(1./Zg2))
    Fw2     = 1.+1j*np.sqrt(np.pi)*w2*sp.erfcx(-1j*w2)
    Q2      = Rp2 + (1. - Rp2)*Fw2
    
    p1      = np.exp(1j*kfreq*(R1+R3))/np.sqrt((R3**3)*R1*(R1+R3))
    p2      = np.exp(1j*kfreq*(R1+R4))/np.sqrt((R4**3)*R1*(R1+R4))
    p3      = np.exp(1j*kfreq*(R2+R3))/np.sqrt((R3**3)*R2*(R2+R3))
    p4      = np.exp(1j*kfreq*(R2+R4))/np.sqrt((R4**3)*R2*(R2+R4))
    
    Int       = (p1 + Q2*p2 + Q1*p3 + Q1*Q2*p4)
    I       = np.sum(Int)*dz
    p_t     = np.sqrt(8.*np.pi*kfreq)*d2*np.exp(-1j*np.pi/4.)/(16.*np.pi**2)*I
    p_t     = np.asscalar(p_t) # put the array back into scalar type for the main script
    R_dir   = np.sqrt(abs(h_s-h_r)**2 + (d1 + d2)**2)
    p_f     = np.exp(1j*kfreq*R_dir)/(4.*np.pi*R_dir)
    return p_f,p_t


def analytic_solution_modes(dx, p, nx, ny, x, y, c, t, n):
    """
     Analytic solution for acoustic modes in 2D rectangular domains.

    :param  dx      spatial step after discretization, scalar (m).
    :param  p       numerical acoustic pressure used for the sizes, 2D-array (Pa).
    :param  nx      mode number (int), scalar.
    :param  ny      mode number (int), scalar.
    :param  x       discrete length sequence of the domain side along the x-axis, scalar (m).
    :param  y       discrete length sequence of the domain side along the y-axis, scalar (m).
    :param  c       sound speed, scalar (m.s-1).
    :param  t       discretized time sequence, 1d array (s).
    :param  n       discrete iteration inside the for loop, scalar.

    :param      xv      vector used for the dimension of the x-axis of the function I_exact, vector.
    :param      yv      vector used for the dimension of the y-axis of the function I_exact, vector.
    :param      NX,NY,Lx_disc,Ly_disc resize the solution according to the exact shape of the inner domain.
    :param      p_exact exact solution for the acoustic pressure modes, 2D-array (Pa).

    :return     the exact solution for the acoustic pressure modes, 2D-array (Pa).
    """
    xv = x[:, np.newaxis] - (dx / 2)  # for vectorized function evaluations
    yv = y[np.newaxis, :] - (dx / 2)

    NX = p.shape[0]-2
    NY = p.shape[1]-2

    Lx_disc = NX * dx
    Ly_disc = NY * dx

    if n==5:
        print '++++++++++++++++++++++++++++++++++++++++'
        print 'shape0=%g ; shape1=%g ' % (p.shape[0], p.shape[1])
        print 'Inner computation nodes: NX=%i, NY=%i.' % (NX, NY)
        print 'Inner lengths: LX=%f, LY=%f.' % (Lx_disc, Ly_disc)
        print '++++++++++++++++++++++++++++++++++++++++'

    def O_exact(nx, ny):
        """Omega, for a rigid box"""
        return np.pi * c * np.sqrt((nx / Lx_disc)**2 +
                                   (ny / Ly_disc)**2)

    def I_exact(xxx, yyy, nx, ny):
        """Spatial solution for a rigid box"""
        return np.cos(np.pi * nx * xxx / Lx_disc) *\
               np.cos(np.pi * ny * yyy / Ly_disc)

    return I_exact(xv[1:-1,:],yv[:,1:-1],nx,ny)*np.cos(O_exact(nx, ny)*t[n])


def analytic_solution_modes_fd_helmholtz(dx, nx, ny, Lx, Ly, Nx, Ny, x, y):
    """
    Analytic solution for acoustic modes in 2D rectangular domains.

    :param dx: spatial step (m).
    :type dx: float
    :param nx: mode number following the x direction
    :type nx: int
    :param ny: mode number following the y direction
    :type ny: int
    :param Lx: length of the box followinf the x axis (m)
    :type Lx: float
    :param Ly: length of the box followinf the y axis (m)
    :type Ly: float
    :param Nx: length of the domain in number of node following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of node following the y dir.
    :type Ny: int
    :return: the analytic solution for acoustic modes in 2D rectangular domains.
    :rtype: 2d-array of floats
    """
    xv = x[:, np.newaxis] - (dx / 2)  # for vectorized function evaluations
    yv = y[np.newaxis, :] - (dx / 2)

    return np.cos(np.pi * nx * xv[1:-1, :] / Lx) * \
           np.cos(np.pi * ny * yv[:, 1:-1] / Ly)

def green_sol_to_2D_helmholtz_eq(Nx, Ny, xs_idx, ys_idx, dx, k):
    """
    Grenn solution for a point source pulse (Dirac) in the 2D Helmholtz
    equation.
    :param Nx: length of the domain in number of node following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of node following the y dir.
    :type Ny: int
    :param xs_idx: node location of the point source for the x coordinate
    :type xs_idx: int
    :param ys_idx: node location of the point source for the y coordinate
    :type ys_idx: int
    :param dx: spatial step (m).
    :type dx: float
    :param k: wavenumber (rad.m-1).
    :type k: float
    :return: the pressure map of the solution.
    :rtype: 2D-numpy array of complex numbers
    """
    green_2d = np.zeros((Nx + 1, Ny + 1), dtype=np.complex128)
    for i in range(1, Nx):
        for j in range(1, Ny):
            green_2d[i, j] = 1j / 4. * sp.hankel1(0, k * np.sqrt(
                                    np.abs((0.001 * dx + xs_idx * dx) - i * dx) +
                                    np.abs((0.001 * dx + ys_idx * dx) - j * dx)))
    return green_2d


def testcase1_eigenfunct(Lx, Ly, dx, x, y):
    """
    Pressure solution for the test case 1 found in [sutmann_jcam2007].

    :param Lx: length of the box followinf the x axis (m)
    :type Lx: float
    :param Ly: length of the box followinf the y axis (m)
    :type Ly: float
    :param Nx: length of the domain in number of node following the x dir.
    :type Nx: int
    :param Ny: length of the domain in number of node following the y dir.
    :type Ny: int
    :param dx: spatial step (m).
    :type dx: float
    :return:
    :rtype:
    """

    xv = x[:, np.newaxis] - (dx / 2)  # for vectorized function evaluations
    yv = y[np.newaxis, :] - (dx / 2)
    return np.sin(np.pi * xv[1:-1, :] / Lx) * np.sin(np.pi * yv[:, 1:-1] / Ly)


def analytic_solution_scattered_pressure(rho, c, f, a, r_rcp, Z, alpha):
    """
     Analytic solution for the pressure scattered by a circular obstacle.

    :param  rho     air density, float (kg.m-3).
    :param  c       sound speed, float (m.s-1).
    :param  f       frequency, float (Hz).
    :param  a       radius of the circular obstacle, float (m).
    :param  r_rcp   distance betwen the center and the reciever, float (m).
    :param  Z       surface impedance of the circular obstacle, flaot ().
    :param  alpha   angle btewen the propagation axis and the receiver-center axis, float (rad).

    :param      omega   angular frequency, float (rad.s-1).
    :param      k       wave number, float (rad.m-1).
    :param      ps      scattered pressure, float (Pa).

    :return     the absolute value of the scattered pressure.
    """
    omega = 2 * np.pi * f
    k = omega / c
    ka1 = k * a
    ps = 0
    for j, order in enumerate(range(8)):
        if order == 0:
            dirac1 = 1
        else:
            dirac1 = 0

        ps = ps - (2 - dirac1) * 1j ** order * \
                  (1j * ((order / ka1) * sp.jn(order, ka1) - sp.jn(order + 1, ka1)) +
                   (rho * c * (Z ** -1)) * sp.jn(order, ka1)) * \
                  (1j * ((order / ka1) * sp.hankel1(order, ka1) - sp.hankel1(order + 1, ka1)) +
                   (rho * c * (Z ** -1)) * sp.hankel1(order, ka1)) ** -1 * \
                  sp.hankel1(order, k * r_rcp) * np.cos(order * alpha)
    return np.abs(ps)


def analytic_solution_slope(freq, h_s, l1max, l2max, hmax, alpha, dx):
    """
      Calculation of the analytic pressure for each point of both
            discretized parts, i.e. horizontal and upward-slopping,
            based on [rasmussen_jsv1985].

    :param      freq    frequency of the source signal, float (Hz).
    :param      h_s     source heights, float (m).
    :param      l1max  length of the horizontal part, float (m).
    :param      l2max  length of the upward-sloping part, float (m).
    :param      hmax    height of the the two parts, float (m).
    :param      alpha   angle of the upward-sloping part, float (rad).
    :param      dx      spatial step, float (m).

    :param  c       sound speed (m.s-1), scalar.
    :param  lamb    wave length (m), scalar.
    :param  k       wave number (rad.m-1), scalar.
    :param  a       calculation parameter used in the v function, scalar.
    :param  alpha   angle of the upward-sloping part (rad), scalar.
    :param  gamma   angle between the two parts following a counterclockwise
                        direction (rad), scalar.
    :param  nu      angular variable used in the v function (rad), scalar.

    :param  dim1    number of points for the horizontal part, scalar.
    :param  dim2    number of points for the upward-sloping part, scalar.
    :param  dimh    number of points for the height of both parts, scalar.

    :param  dx1     spatial step follwing the length of the horizontal part, scalar.
    :param  dx2     spatial step follwing the length of the upward-sloping part, scalar.
    :param  dxh     spatial step follwing the height of both parts, scalar.

    :param  l1      discretized length of the horizontal part, array.
    :param  l2      discretized length of the upward-sloping part, array.
    :param  h       discretized heights of both parts, array.

    :param  r1,r2,r3,r4 the four paths, scalar.
    :param  p0      pressure field of the horizontal part, complex.
    :param  pdiff0  diffracted pressure field of the horizontal part, complex.
    :param  p1      pressure field of the upward-slopping part, complex.
    :param  pdiff1  diffracted pressure field of the upward-slopping part, complex.
    :param  res01   absolute value of the incident pressure for the flat part, float (Pa).
    :param  res02   absolute value of the diffracted pressure for the flat part, float (Pa).
    :param  res03   absolute value of the total pressure: incident+diffracted for the flat part, float (Pa).
    :param  res11   absolute value of the incident pressure for the upward-slopping part, float (Pa).
    :param  res12   absolute value of the diffracted pressure for the upward-slopping part, float (Pa).
    :param  res13   absolute value of the total pressure: incident+diffracted for the upward-slopping part, float (Pa).
    :return The absolute values of the pressure fields for each point of
            each domain in the 2D-arrays res01 to res13.
    """
    c = 340.
    lamb = c / freq
    k = 2. * np.pi / lamb
    alpha = np.pi * alpha / 180.
    gamma = np.arctan(h_s / l1max)
    nu = 2. - (np.pi + alpha) / np.pi

    l1 = np.arange(0, l1max + dx, dx)
    l2 = np.arange(0, l2max + dx, dx)
    h = np.arange(dx, hmax + dx, dx)
    # ==============================================================================
    #    # ------- Horizontal part (dim1*dimh), before corner ------- #
    # ==============================================================================
    res01 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    res02 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    res03 = np.zeros((len(l1), len(h)), dtype=np.complex128)
    for mm in range(1, len(l1)):
        #        print dim1-mm
        beta0 = np.arctan(h / (l1max - l1[mm]))
        beta = np.pi - alpha - beta0

        r1 = np.sqrt(l1[mm] ** 2 + (h - h_s) ** 2)
        r2 = np.sqrt(l1[mm] ** 2 + (h + h_s) ** 2)

        p0 = (1. / 4) * np.sqrt(2. / (np.pi * k * r1)) * (np.exp(1j * k * r1)) \
             + (1. / 4) * np.sqrt(2. / (np.pi * k * r2)) * (np.exp(1j * k * r2))

        # ------- Diffracted part ------- #
        r0diff = np.sqrt(h_s ** 2 + l1max ** 2)
        r1diff = np.sqrt((l1max - l1[mm]) ** 2 + h ** 2)
        gr1diff = r0diff + r1diff
        pdiff0 = (1. / 4) * np.sqrt(2. / (np.pi * k * gr1diff)) * np.exp(1j * k * (gr1diff)) \
                 * (v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta + gamma, nu, len(h), k)
                    + v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta - gamma, nu, len(h), k))

        res01[mm - 1, :] = np.abs(p0)
        res02[mm - 1, :] = np.abs(pdiff0)
        res03[mm - 1, :] = np.abs(p0 + pdiff0)

    # ==============================================================================
    #      # ------- Upward-sloping part (dim2*dimh) after corner ------- #
    # ==============================================================================
    res11 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    res12 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    res13 = np.zeros((len(l2), len(h)), dtype=np.complex128)
    for mm in range(1, len(l2)):
        #        print dim2-mm
        beta = np.asarray([np.pi / 2 for i in range(len(h))])
        if mm > 2:
            beta = np.arctan(h / l2[mm])

        r0diff = np.sqrt(h_s ** 2 + l1max ** 2)
        r1diff = np.sqrt(l2[mm] ** 2 + h ** 2)
        gr1diff = r0diff + r1diff
        r1 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha + gamma + beta))
        r2 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha - gamma + beta))
        r3 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha + gamma - beta))
        r4 = np.sqrt(l1max ** 2 + r1diff ** 2 + 2. * l1max * r1diff * np.cos(alpha - gamma - beta))
        h4 = l2[mm] * np.tan(alpha - gamma)
        h3 = l2[mm] * np.tan(alpha + gamma)
        p1 = np.zeros((len(h)), dtype=np.complex128)
        for ll in range(len(h)):
            p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r1[ll])) * (np.exp(1j * k * r1[ll])) \
                      + (1. / 4) * np.sqrt(2. / (np.pi * k * r2[ll])) * (np.exp(1j * k * r2[ll]))
            if h[ll] <= h3:
                p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r3[ll])) * (np.exp(1j * k * r3[ll]))
            if h[ll] <= h4:
                p1[ll] += (1. / 4) * np.sqrt(2. / (np.pi * k * r4[ll])) * (np.exp(1j * k * r4[ll]))

        # ------- Diffracted part ------- #
        pdiff1 = (1. / 4) * np.sqrt(2. / (np.pi * k * gr1diff)) * np.exp(1j * k * (gr1diff)) \
                 * (v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta + gamma, nu, len(h), k)
                    + v(r0diff * r1diff / gr1diff, 1., np.pi - alpha - beta - gamma, nu, len(h), k))

        res11[mm - 1, :] = np.abs(p1)
        res12[mm - 1, :] = np.abs(pdiff1)
        res13[mm - 1, :] = np.abs(p1 + pdiff1)
    return res03, res13


def fp(x):
    """Function used in **v(a, b, th, nu, dimh, k)** for **analytic_solution_slope()**
    :param x: real number
    :type x: list
    :return: fp value
    :rtype: list
    """
    rx = np.sqrt(x * 2 / np.pi)
    s_fresnel, c_fresnel = sp.fresnel(rx)
    return - 2 * 1j * np.sqrt(x) * np.exp(-1j * x) * np.sqrt(np.pi / 2.) \
           * (.5 - c_fresnel + 1j * (.5 - s_fresnel))


def v(a, b, th, nu, dimh, k):
    """Function used in **analytic_solution_slope()**
    :param a:
    :type a:
    :param b:
    :type b:
    :param th:
    :type th:
    :param nu:
    :type nu:
    :param dimh:
    :type dimh:
    :param k:
    :type k:
    :return:
    :rtype:
    """
    #    real, b
    #    real, dimension(dimh) :: a, th
    #    real, dimension(dimh) :: vp, vm, v
    #    integer, dimension(dimh) :: nplus, nmoins
    vp = np.zeros((dimh), dtype=np.complex128)
    vm = np.zeros((dimh), dtype=np.complex128)
    nmoins = np.zeros((dimh), dtype=np.float64)
    nplus = np.ones((dimh), dtype=np.float64)
    for th_idx, th_val in enumerate(th):
        if (th_val - nu * np.pi + np.pi) <= 0:
            nplus[th_idx] = 0
        nmoins[th_idx] = 1
        if (th_val - nu * np.pi - np.pi) <= 0:
            nmoins[th_idx] = 0
        if (th_val + nu * np.pi - np.pi) <= 0:
            nmoins[th_idx] = -1
        vp[th_idx] = -  (np.exp(1j * np.pi / 4)) * (2 * np.pi * k * a[th_idx] * b) ** (-.5) / (2 * nu) * \
                     (np.tan((np.pi + th_val) / (2 * nu))) ** (-1) * \
                     fp(2 * k * a[th_idx] * (np.cos((2 * nplus[th_idx] * nu * np.pi - th_val) / 2)) ** 2)
        vm[th_idx] = -  (np.exp(1j * np.pi / 4)) * (2 * np.pi * k * a[th_idx] * b) ** (-.5) / (2 * nu) * \
                     (np.tan((np.pi - th_val) / (2 * nu))) ** (-1) * \
                     fp(2 * k * a[th_idx] * (np.cos((2 * nmoins[th_idx] * nu * np.pi - th_val) / 2)) ** 2)
    return vp + vm