# -*- coding: utf-8 -*-
##
# \file     pe_absorbing_layer.py
# \title    Definition of an absorbing layer for the parabolic equation.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 20 Nov.
##
import numpy as np


def abs_lay_top(p_ij, dy, Ny, y_start_abs):
    """
    Absorbing layer for the parabolic equation defined following the
    vertical direction - i.e. along the y axis.
    For more details see **[chevret_phd1994, Eq.(4.38), p.59]**.

    :param p_ij: pressure at the discrete location i,j ~ (x, y) (Pa).
    :type p_ij: 2D numpy arrays of complexes
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :param y_start_abs: y coordinate of the layer starting position (m).
    :type y_start_abs: float
    :return: the pressure array inclunding the absorption of the layer (Pa).
    :rtype: 2D numpy arrays of complexes
    """
    a_empirical = 4.5
    b_empirical = int(round(1.4/dy))
    abs_lay_top = np.ones((Ny + 1), dtype=np.float64)
    for j in range(int(y_start_abs / dy) + 1, Ny + 1):
        abs_lay_top[j] = np.exp(-((j - int(y_start_abs / dy)) /
                                (a_empirical * (Ny + b_empirical - j)))**2)
        p_ij[j] = p_ij[j] * abs_lay_top[j]

    plot_abs_profil = False
    if plot_abs_profil:
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.plot(range(int(y_start_abs / dy) + 1, Ny + 1),
                 abs_lay_top[int(y_start_abs / dy) + 1: Ny + 1])
        plt.show()
    return p_ij


def abs_lay_bottom_top(p_ij, dy, Ny, y_start_abs):
    """
    Absorbing layer for the ground (low part) of the parabolic equation, in
    order to simulate free field propagation.
    The absorbing layer is defined following the vertical direction -
    i.e. along the y axis.
    For more details see **[chevret_phd1994, Eq.(4.38), p.59]**.

    :param p_ij: pressure at the discrete location i,j ~ (x, y) (Pa).
    :type p_ij: 2D numpy arrays of complexes
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :param y_start_abs: y coordinate of the layer starting position (m).
    :type y_start_abs: float
    :return: the pressure array inclunding the absorption of the layer (Pa).
    :rtype: 2D numpy arrays of complexes
    """
    a_empirical = 4.5
    b_empirical = int(round(1.4/dy))
    abs_lay_top = np.ones((Ny + 1), dtype=np.float64)
    abs_lay_total = np.ones((Ny + 1), dtype=np.float64)
    abs_lay_bottom = np.ones((Ny + 1), dtype=np.float64)
    for j in range(int(y_start_abs / dy) + 1, Ny + 1):
        abs_lay_top[j] = np.exp(-((j - int(y_start_abs / dy)) /
                                (a_empirical * (Ny + b_empirical - j)))**2)
        abs_lay_bottom[Ny - j - 1] = np.exp(-((j + 1 - int(y_start_abs / dy)) /
                                              (a_empirical *
                                               (- Ny - b_empirical + j))) ** 2)

    for j in range(Ny + 1):
        abs_lay_total[j] = abs_lay_bottom[j] * abs_lay_top[j]
        p_ij[j] = p_ij[j] * abs_lay_total[j]

    plot_abs_profil = False
    if plot_abs_profil:
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.plot(range(int(y_start_abs / dy) + 1, Ny + 1),
                 abs_lay_top[int(y_start_abs / dy) + 1: Ny + 1])
        plt.figure(2)
        plt.plot(range(Ny - int(y_start_abs / dy)),
                 abs_lay_bottom[: Ny - int(y_start_abs / dy)])
        plt.show()

    return p_ij


def abs_lay_top_1(k, dy, Ny, y_start_abs):
    """
    Absorbing layer for the parabolic equation defined following the
    vertical direction - It has to be applied on the wavenumber directly.

    :param p_ij: pressure at the discrete location i,j ~ (x, y) (Pa).
    :type p_ij: 2D numpy arrays of complexes
    :param dy: spatial step for the y directions (m).
    :type dy: float
    :param Ny: length of the domain in number of nodes following the y dir.
    :type Ny: int
    :param y_start_abs: y coordinate of the layer starting position (m).
    :type y_start_abs: float
    :return: the wavenumber.
    :rtype: 2D numpy arrays of complexes
    """
    abs_lay_top = np.ones((Ny + 1), dtype=np.float64)
    A = np.ones((Ny + 1), dtype=np.float64)
    for j in range(int(y_start_abs / dy) + 1, Ny + 1):
        A[Ny + 1 - j + int(y_start_abs / dy)] = np.exp(-4. *
                      ((j * dy) - y_start_abs) / ((Ny + 1) * dy - y_start_abs))
        abs_lay_top[Ny + 1 - j + int(y_start_abs / dy)] = A[j] * \
                              ((j * dy - y_start_abs) /
                               ((Ny + 1) * dy - y_start_abs)) ** 2
        k[j] = k[j] * abs_lay_top[j]
    return k
