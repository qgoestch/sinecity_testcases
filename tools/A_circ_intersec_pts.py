# -*- coding: utf-8 -*-
##
# \file     A_circ_intersec_pts.py
# \title    Fetch the intersection points between the continuous circle and the grid.
#           A_circ is the first step over four for the definition of a circular obstacle
#           in the FVDT, FDTD and TLM methods.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 07 Sep.
##

import numpy as np
import math

def intesection_points(ray, dx, ax, by):
    """
          Step 1/4, fetch the instersection points between the continuous circle and the grid.
    The interescetion points are fetched quadrant per quadrant, starting on the upper-right side.

    :param  ray     radius of the circle, float (m).
    :param  dx      spatial step set, float (m).
    :param  ax      x coordinate of the circle center, float (m).
    :param  by      y coordinate of the circle center, float (m).

    :param      x_inter x coordinates of the interesection points, list of floats (m).
    :param      y_inter y coordinates of the interesection points, list of floats (m).
    :param      x_integ intermediate x coordinates to fetch the inters. pts, list of floats (m).
    :param      y_integ intermediate y coordinates to fetch the inters. pts, list of floats (m).
    :param      x_ins   x coordinates of the nearest grid points inside the circle, list of floats (m).
    :param      y_ins   y coordinates of the nearest grid points inside the circle, list of floats (m).
    :param      x_out   x coordinates of the nearest grid points outside the circle, list of floats (m).
    :param      y_out   y coordinates of the nearest grid points outside the circle, list of floats (m).

    :return     the instersection points x_inter, y_inter.
    """
    xs = int(np.around(ray / dx * np.cos(0))) * dx + ax  # starting points
    ys = int(np.around(ray / dx * np.sin(0))) * dx + by
    alpha = []
    x_inter = [];
    y_inter = []
    x_integ = [];
    y_integ = []
    x_ins = [];
    y_ins = []
    x_out = [];
    y_out = []
    # ==============================================================================
    # Initial shift: dx/2 as xs & ys are cell CENTER
    # ==============================================================================
    yi = ys + dx / 2  # closest y-point
    x_new = np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection
    alpha_x = np.tan((yi - by) / (x_new - ax))  # new angle on the cirlcle

    xi = xs - dx / 2  # closest x-point
    y_new = np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection
    alpha_y = np.tan((y_new - by) / (xi - ax))  # new angle on the cirlcle

    k = 0
    if alpha_x < alpha_y:
        k += 1
        x_inter.append(x_new)
        y_inter.append(yi)
        x_ins.append(xi)
        y_ins.append(yi)
        x_out.append(xs + dx / 2)
        y_out.append(yi)
        x_integ.append(xs + dx / 2)
        y_integ.append(yi)
        alpha.append(alpha_x)

        j = 0
        while True:
            yi = y_integ[j] + dx  # closest y-point
            x_new = np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection
            alpha_x = np.tan((yi - by) / (x_new - ax))  # new angle on the cirlcle

            xi = x_integ[j] - dx  # closest x-point
            y_new = np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection
            alpha_y = np.tan((y_new - by) / (xi - ax))  # new angle on the cirlcle

            if alpha_x < alpha_y:
                x_inter.append(x_new)
                y_inter.append(yi)
                x_ins.append(xi)
                y_ins.append(yi)
                x_out.append(xi + dx)
                y_out.append(yi)
                x_integ.append(x_integ[j])
                y_integ.append(yi)
                alpha.append(alpha_x)

            if alpha_y < alpha_x:
                x_inter.append(xi)
                y_inter.append(y_new)
                x_ins.append(xi)
                y_ins.append(yi - dx)
                x_out.append(xi)
                y_out.append(yi)
                x_integ.append(xi)
                y_integ.append(y_integ[j])
                alpha.append(alpha_y)

            j += 1
            if j == 2: break

    if k == 0:
        import sys
        sys.exit('Pb, exit, h=%f' % dx)

    i = j - 1
    count2 = 0
    count3 = 0
    count4 = 0
    while ((x_inter[i + 1] != x_inter[0]) or (y_inter[i + 1] != y_inter[0])):
        i += 1
        #       1st quadrant
        if ((y_integ[i] - by) > 0 and (x_integ[i] - ax) > 0):

            yi = y_integ[i] + dx  # closest y-point
            x_new = np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection

            xi = x_integ[i] - dx  # closest x-point
            y_new = np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection

            if math.isnan(y_new) or ((yi - by) / (x_new - ax)) < ((y_new - by) / (xi - ax)):
                x_inter.append(x_new)
                y_inter.append(yi)
                x_ins.append(xi)
                y_ins.append(yi)
                x_out.append(xi + dx)
                y_out.append(yi)
                x_integ.append(x_integ[i])
                y_integ.append(yi)

            if math.isnan(x_new) or ((yi - by) / (x_new - ax)) > ((y_new - by) / (xi - ax)):
                x_inter.append(xi)
                y_inter.append(y_new)
                x_ins.append(xi)
                y_ins.append(yi - dx)
                x_out.append(xi)
                y_out.append(yi)
                x_integ.append(xi)
                y_integ.append(y_integ[i])

        #       2nd quadrant
        if ((y_integ[i] - by) > 0) and ((x_integ[i] - ax) < 0):
            count2 += 1
            if count2 == 1:    y_integ[i] = y_integ[i] + dx

            yi = y_integ[i] - dx  # closest y-point
            x_new = -np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection

            xi = x_integ[i] - dx  # closest x-point
            y_new = np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection

            if math.isnan(y_new) or ((yi - by) / (x_new - ax)) < ((y_new - by) / (xi - ax)):
                x_inter.append(x_new)
                y_inter.append(yi)
                x_ins.append(xi + dx)
                y_ins.append(yi)
                x_out.append(xi)
                y_out.append(yi)
                x_integ.append(x_integ[i])
                y_integ.append(yi)

            if math.isnan(x_new) or ((yi - by) / (x_new - ax)) > ((y_new - by) / (xi - ax)):
                x_inter.append(xi)
                y_inter.append(y_new)
                x_ins.append(xi)
                y_ins.append(yi)
                x_out.append(xi)
                y_out.append(yi + dx)
                x_integ.append(xi)
                y_integ.append(y_integ[i])

        #       3rd quadrant
        if ((y_integ[i] - by) < 0) and ((x_integ[i] - ax) < 0):
            count3 += 1
            if count3 == 1:    x_integ[i] = x_integ[i] - dx
            #            print '3rd quadrant'
            #            break
            yi = y_integ[i] - dx  # closest y-point
            x_new = -np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection

            xi = x_integ[i] + dx  # closest x-point
            y_new = -np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection

            if math.isnan(y_new) or ((yi - by) / (x_new - ax)) < ((y_new - by) / (xi - ax)):
                x_inter.append(x_new)
                y_inter.append(yi)
                x_ins.append(xi)
                y_ins.append(yi)
                x_out.append(xi - dx)
                y_out.append(yi)
                x_integ.append(x_integ[i])
                y_integ.append(yi)

            if math.isnan(x_new) or ((yi - by) / (x_new - ax)) > ((y_new - by) / (xi - ax)):
                x_inter.append(xi)
                y_inter.append(y_new)
                x_ins.append(xi)
                y_ins.append(yi + dx)
                x_out.append(xi)
                y_out.append(yi)
                x_integ.append(xi)
                y_integ.append(y_integ[i])


        #       4th quadrant
        if ((y_integ[i] - by) < 0) and ((x_integ[i] - ax) > 0):
            count4 += 1
            if count4 == 1:    y_integ[i] = y_integ[i] - dx
            #            print '4th quadrant'
            yi = y_integ[i] + dx  # closest y-point
            x_new = np.sqrt(ray ** 2 - (yi - by) ** 2) + ax  # new x-coordinate = intersection

            xi = x_integ[i] + dx  # closest x-point
            y_new = -np.sqrt(ray ** 2 - (xi - ax) ** 2) + by  # new y-coordinate = intersection

            if math.isnan(y_new) or ((yi - by) / (x_new - ax)) < ((y_new - by) / (xi - ax)):
                x_inter.append(x_new)
                y_inter.append(yi)
                x_ins.append(xi - dx)
                y_ins.append(yi)
                x_out.append(xi)
                y_out.append(yi)
                x_integ.append(x_integ[i])
                y_integ.append(yi)

            if math.isnan(x_new) or ((yi - by) / (x_new - ax)) > ((y_new - by) / (xi - ax)):
                x_inter.append(xi)
                y_inter.append(y_new)
                x_ins.append(xi)
                y_ins.append(yi)
                x_out.append(xi)
                y_out.append(yi - dx)
                x_integ.append(xi)
                y_integ.append(y_integ[i])

    # ==============================================================================
    # ------------ Plot GRID vs. CIRCLE with the intersection points ---------------
    # ==============================================================================
    plot_intersections=False
    if plot_intersections:
        phi = np.linspace(0, 2 * np.pi, 10000)
        x_circ = np.zeros((len(phi)), dtype=np.float64)
        x_circ = ray * np.cos(phi) + ax
        y_circ = np.zeros((len(phi)), dtype=np.float64)
        y_circ = ray * np.sin(phi) + by

        import matplotlib.pyplot as plt
        fig = plt.figure('GRID vs. CIRCLE (1)')
        ax1  = fig.add_subplot(111)
        ax1.set_aspect('equal')

        plt.plot(ax,by,'ko', markersize=10, markeredgewidth=4, markeredgecolor='k', markerfacecolor='None')
        plt.plot(x_circ,y_circ,'c-')
        plt.plot(x_inter,y_inter,'ro--')

        plt.plot(x_ins, y_ins, 'gs--')
        plt.plot(x_out, y_out, 'md--')

        ax1.xaxis.set_ticks(np.arange(0,2*ax+dx,dx))
        ax1.yaxis.set_ticks(np.arange(0,2*by+dx,dx))
        plt.grid(True, which='major', color='b', linestyle='--',lw=1.)

        plt.minorticks_on()
        ax1.set_xticks(np.arange(dx/2,2*ax+dx+dx/2,dx), minor=True)
        ax1.set_yticks(np.arange(dx/2,2*by+dx+dx/2,dx), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-',lw=3.)
        radius = ray
        plt.xlim(ax- 1.2*radius, ax+ 1.4*radius)
        plt.ylim(by- 1.2*radius, by+ 1.4*radius)
        plt.show()

    return x_inter, y_inter, x_integ, y_integ, x_ins, y_ins, x_out, y_out