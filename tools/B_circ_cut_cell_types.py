# -*- coding: utf-8 -*-
##
# \file     B_circ_cut_cell_types.py
# \title    Get the parameters of the cut cells.
#           SWEEP ACROSS ALL INTERSECTION POINTS TO CALCULATE
#           THE NEW VOLUMES, CUT SURFACES, SIDE SURFACES, COORDINATES OF THE CENTERS OF MASS.
#           This parameters are required for FVTD calculations only, but is mandatory to get the corner
#           of the FDTD circle.
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 07 Sep.
##

import numpy as np
import A_circ_intersec_pts as pts

def cut_cell_types(radius, h, ax, by):
    """
          Step 2/4, calculation of the cut-cells parameters for the FVTD update.

    :param  radius  radius of the circle, float (m).
    :param  h       spatial step set, float (m).
    :param  ax      x coordinate of the circle center, float (m).
    :param  by      y coordinate of the circle center, float (m).

    :param      x_halcell   x coordinates of the halved cells, list of floats (m).
    :param      x_halcell   y coordinates of the halved cells, list of floats (m).
    :param      x_cutcell   x coordinates of the cut cells, list of floats (m).
    :param      y_cutcell   y coordinates of the cut cells, list of floats (m).
    :param      x_cc   x coordinates of the new center of mass for the cut cells, list of floats (m).
    :param      y_cc   y coordinates of the new center of mass for the cut cells, list of floats (m).

    :return     among intermediate parameters, it resturns the coordinates of the boundary cells (halcell or cutcell).
    """
    x_inter, y_inter, x_integ, y_integ, \
    x_ins, y_ins, x_out, y_out = pts.intesection_points(radius, h, ax, by)
    # ==============================================================================
    #     intesection_points(radius:ray, spatial_step:h, x_center:ax, y_center:by)
    # ==============================================================================
    x_case1 = [];
    y_case1 = [];
    x_case9 = [];
    y_case9 = [];
    x_case2 = [];
    y_case2 = [];
    x_case10 = [];
    y_case10 = [];
    x_case3 = [];
    y_case3 = [];
    x_case11 = [];
    y_case11 = [];
    x_case4 = [];
    y_case4 = [];
    x_case12 = [];
    y_case12 = [];
    x_case5 = [];
    y_case5 = [];
    x_case13 = [];
    y_case13 = [];
    x_case6 = [];
    y_case6 = [];
    x_case14 = [];
    y_case14 = [];
    x_case7 = [];
    y_case7 = [];
    x_case15 = [];
    y_case15 = [];
    x_case8 = [];
    y_case8 = [];
    x_case16 = [];
    y_case16 = [];

    x_cc = [];
    y_cc = []
    S = [];
    V = [];
    S_1 = [];
    S_2 = [];
    S_3 = [];
    S_4 = [];
    case = []
    # ==============================================================================
    #       SWEEP ACROSS ALL INTERSECTION POINTS TO CALCULATE
    #       FOR ALL CUT-CELLS:
    #       - THE VOLUMES (V) AND CUT-SURFACES (S),
    #       - THE SIDE SURFACES: S_1; S_2; S_3; S_4;
    #       - THE COORDINATES OF THE CENTERS: x_cc, y_cc
    # ==============================================================================
    for i in range(len(x_inter) - 1):
        #       1st quadrant
        if ((y_integ[i] - by) > 0 and (x_integ[i] - ax) > 0):
            if i == len(x_inter) - 1:
                j = 0
            else:
                j = i + 1
            if (x_inter[i] != x_integ[i] and x_inter[j] != x_integ[j] and
                        y_inter[i] == y_integ[i] and y_inter[j] == y_integ[j]):
                #                   CASE 1
                case.append(1)
                x_case1.append(x_integ[i] - h / 2)
                y_case1.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] + h / 2)
                S.append(np.sqrt(h ** 2 + (x_inter[i] - x_inter[i + 1]) ** 2))
                V.append(h ** 2 - (0.5 * h * (x_inter[i] - x_inter[i + 1]) +
                                   h * (x_inter[i + 1] - (x_integ[i + 1] - h))))
                S_1.append(x_integ[i + 1] - x_inter[i + 1])
                S_2.append(h)
                S_3.append(x_integ[i] - x_inter[i])
                S_4.append(0)

            if (x_inter[i] != x_integ[i] and x_inter[j] == x_integ[j] and
                        y_inter[i] == y_integ[i] and y_inter[j] != y_integ[j]):
                #                   CASE 2
                case.append(2)
                x_case2.append(x_integ[i] - h / 2)
                y_case2.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] + h / 2)
                S.append(np.sqrt((x_inter[i] - x_integ[i + 1]) ** 2 +
                                 (y_inter[i + 1] - y_integ[i]) ** 2))
                V.append(h ** 2 - 0.5 * (x_inter[i] - x_integ[i + 1]) *
                         (y_inter[i + 1] - y_integ[i]))
                S_1.append(h)
                S_2.append(h)
                S_3.append(h - (x_inter[i] - x_integ[i + 1]))
                S_4.append(h - (y_inter[i + 1] - y_integ[i]))

            if (x_inter[i] == x_integ[i] and x_inter[j] != x_integ[j] and
                        y_inter[i] != y_integ[i] and y_inter[j] == y_integ[j]):
                #                   CASE 3
                case.append(3)
                x_case3.append(x_integ[i] - h / 2)
                y_case3.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] + h / 2)
                S.append(np.sqrt((x_integ[i + 1] - x_inter[i + 1]) ** 2 +
                                 (y_integ[i + 1] - y_inter[i]) ** 2))
                V.append(0.5 * (x_inter[i] - x_inter[i + 1]) *
                         (y_inter[i + 1] - y_inter[i]))

                S_1.append(x_integ[i + 1] - x_inter[i + 1])
                S_2.append(y_integ[i + 1] - y_inter[i])
                S_3.append(0)
                S_4.append(0)

            if (x_inter[i] == x_integ[i] and x_inter[j] == x_integ[j] and
                        y_inter[i] != y_integ[i] and y_inter[j] != y_integ[j]):
                #                   CASE 4
                case.append(4)
                x_case4.append(x_integ[i] - h / 2)
                y_case4.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] + h / 2)
                S.append(np.sqrt(h ** 2 + (y_inter[i + 1] - y_inter[i]) ** 2))
                V.append(h ** 2 - (0.5 * h * (y_inter[i + 1] - y_inter[i]) +
                                   h * (y_inter[i] - y_integ[i])))

                S_1.append(h)
                S_2.append(h - (y_inter[i] - y_integ[i]))
                S_3.append(0)
                S_4.append(h - (y_inter[i + 1] - y_integ[i]))

                #       2nd quadrant
        if ((y_integ[i] - by) > 0 and (x_integ[i] - ax) < 0):

            if (x_inter[i] == x_integ[i] and x_inter[i + 1] == x_integ[i + 1] and
                        y_inter[i] != y_integ[i] and y_inter[i + 1] != y_integ[i + 1]):
                #                   CASE 5
                case.append(5)
                x_case5.append(x_integ[i] - h / 2)
                y_case5.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] - h / 2)
                S.append(np.sqrt(h ** 2 + (y_inter[i] - y_inter[i + 1]) ** 2))
                V.append(0.5 * h * (y_inter[i] - y_inter[i + 1]) +
                         h * (y_integ[i + 1] - y_inter[i]))

                S_1.append(h)
                S_2.append(y_integ[i] - y_inter[i])
                S_3.append(0)
                S_4.append(y_integ[i] - y_inter[i + 1])

            if (x_inter[i] == x_integ[i] and x_inter[i + 1] != x_integ[i + 1] and
                        y_inter[i] != y_integ[i] and y_inter[i + 1] == y_integ[i + 1]):
                #                   CASE 6
                case.append(6)
                x_case6.append(x_integ[i] - h / 2)
                y_case6.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] - h / 2)
                S.append(np.sqrt((x_inter[i] - x_inter[i + 1]) ** 2 +
                                 (y_inter[i] - y_inter[i + 1]) ** 2))
                V.append(h ** 2 - (0.5 * (x_inter[i] - x_inter[i + 1]) *
                                   (y_inter[i] - y_inter[i + 1])))

                S_1.append(h)
                S_2.append(h - (y_inter[i] - y_inter[i + 1]))
                S_3.append(h - (x_inter[i] - x_inter[i + 1]))
                S_4.append(h)

            if (x_inter[i] != x_integ[i] and x_inter[i + 1] == x_integ[i + 1] and
                        y_inter[i] == y_integ[i] and y_inter[i + 1] != y_integ[i + 1]):
                #                   CASE 7
                case.append(7)
                x_case7.append(x_integ[i] - h / 2)
                y_case7.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] - h / 2)
                S.append(np.sqrt((x_inter[i] - x_integ[i + 1]) ** 2 +
                                 (y_integ[i + 1] - y_inter[i + 1]) ** 2))
                V.append(0.5 * (x_inter[i] - x_integ[i + 1]) *
                         (y_integ[i + 1] - y_inter[i + 1]))

                S_1.append(x_inter[i] - x_inter[i + 1])
                S_2.append(0)
                S_3.append(0)
                S_4.append(y_integ[i] - y_inter[i + 1])

            if (x_inter[i] != x_integ[i] and x_inter[i + 1] != x_integ[i + 1] and
                        y_inter[i] == y_integ[i] and y_inter[i + 1] == y_integ[i + 1]):
                #                   CASE 8
                case.append(8)
                x_case8.append(x_integ[i] - h / 2)
                y_case8.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] - h / 2)
                y_cc.append(y_integ[i] - h / 2)
                S.append(np.sqrt(h ** 2 + (x_inter[i] - x_inter[i + 1]) ** 2))
                V.append(h ** 2 - (0.5 * h * (x_inter[i] - x_inter[i + 1]) +
                                   h * (x_integ[i] - x_inter[i])))

                S_1.append(h - (x_integ[i] - x_inter[i]))
                S_2.append(0)
                if (x_integ[i + 1] - x_inter[i + 1]) > 0:
                    S_3.append(h - (x_integ[i + 1] - x_inter[i + 1]))
                else:
                    S_3.append(x_inter[i + 1] - x_integ[i + 1])
                S_4.append(h)

                #       3rd quadrant
        if ((y_integ[i] - by) < 0 and (x_integ[i] - ax) < 0):
            if (x_inter[i] != x_integ[i] and x_inter[i + 1] != x_integ[i + 1] and
                        y_inter[i] == y_integ[i] and y_inter[i + 1] == y_integ[i + 1]):
                #                   CASE 9
                case.append(9)
                x_case9.append(x_integ[i] + h / 2)
                y_case9.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] - h / 2)
                V.append(0.5 * h * (x_inter[i + 1] - x_inter[i]) +
                         h * (x_inter[i] - x_integ[i]))
                S.append(np.sqrt((x_inter[i + 1] - x_inter[i]) ** 2 + h ** 2))

                S_1.append(x_inter[i] - x_integ[i])
                S_2.append(0)
                S_3.append(x_inter[i + 1] - x_integ[i + 1])
                S_4.append(h)

            if (x_inter[i] != x_integ[i] and x_inter[i + 1] == x_integ[i + 1] and
                        y_inter[i] == y_integ[i] and y_inter[i + 1] != y_integ[i + 1]):
                #                   CASE 10
                case.append(10)
                x_case10.append(x_integ[i] + h / 2)
                y_case10.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] - h / 2)
                V.append(h ** 2 - (0.5 * (x_inter[i + 1] - x_inter[i]) *
                                   (y_inter[i] - y_inter[i + 1])))
                S.append(np.sqrt((x_integ[i + 1] - x_inter[i]) ** 2 +
                                 (y_inter[i] - y_inter[i + 1]) ** 2))

                S_1.append(h - (x_inter[i + 1] - x_inter[i]))
                S_2.append(h - (y_inter[i] - y_inter[i + 1]))
                S_3.append(h)
                S_4.append(h)

            if (x_inter[i] == x_integ[i] and x_inter[i + 1] != x_integ[i + 1] and
                        y_inter[i] != y_integ[i] and y_inter[i + 1] == y_integ[i + 1]):
                #                   CASE 11
                case.append(11)
                x_case11.append(x_integ[i] + h / 2)
                y_case11.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] - h / 2)
                V.append((0.5 * (x_inter[i + 1] - x_inter[i]) *
                          (y_inter[i] - y_inter[i + 1])))
                S.append(np.sqrt((x_inter[i + 1] - x_inter[i]) ** 2 +
                                 (y_inter[i] - y_inter[i + 1]) ** 2))
                S_1.append(0)
                S_2.append(0)
                S_3.append(x_inter[i + 1] - x_inter[i])
                S_4.append(y_inter[i] - y_inter[i + 1])

            if (x_inter[i] == x_integ[i] and x_inter[i + 1] == x_integ[i + 1] and
                        y_inter[i] != y_integ[i] and y_inter[i + 1] != y_integ[i + 1]):
                #                   CASE 12
                case.append(12)
                x_case12.append(x_integ[i] + h / 2)
                y_case12.append(y_integ[i] - h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] - h / 2)
                V.append(h ** 2 - (0.5 * h * (y_inter[i] - y_inter[i + 1]) +
                                   h * (y_integ[i] - y_inter[i])))
                S.append(np.sqrt(h ** 2 + (y_inter[i] - y_inter[i + 1]) ** 2))

                S_1.append(0)
                if (y_integ[i + 1] - y_inter[i + 1]) > 0:
                    S_2.append(h - (y_integ[i + 1] - y_inter[i + 1]))
                else:
                    S_2.append(y_inter[i + 1] - y_integ[i + 1])
                S_3.append(h)
                S_4.append(h - (y_integ[i] - y_inter[i]))

                #       4th quadrant
        if ((y_integ[i] - by) < 0 and (x_integ[i] - ax) > 0):
            if i == len(x_inter) - 1:
                j = 0
            else:
                j = i + 1

            if (x_inter[i] == x_integ[i] and x_inter[j] == x_integ[j] and
                        y_inter[i] != y_integ[i] and y_inter[j] != y_integ[j]):
                #                   CASE 13
                case.append(13)
                x_case13.append(x_integ[i] + h / 2)
                y_case13.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] + h / 2)
                V.append(0.5 * h * (y_inter[i + 1] - y_inter[i]) +
                         h * (y_inter[i] - y_integ[i]))
                S.append(np.sqrt((x_inter[i + 1] - x_inter[i]) ** 2 +
                                 (y_inter[i + 1] - y_integ[i]) ** 2))

                S_1.append(0)
                S_2.append(y_inter[i + 1] - y_integ[i + 1])
                S_3.append(h)
                S_4.append(y_inter[i] - y_integ[i])

            if (x_inter[i] == x_integ[i] and x_inter[j] != x_integ[j] and
                        y_inter[i] != y_integ[i] and y_inter[j] == y_integ[j]):
                #                   CASE 14
                case.append(14)
                x_case14.append(x_integ[i] + h / 2)
                y_case14.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] + h / 2)
                V.append(h ** 2 - (0.5 * (x_inter[i + 1] - x_inter[i]) *
                                   (y_inter[i + 1] - y_inter[i])))
                S.append(np.sqrt((x_inter[i + 1] - x_inter[i]) ** 2 +
                                 (y_inter[i + 1] - y_inter[i]) ** 2))

                S_1.append(h - (x_inter[i + 1] - x_inter[i]))
                S_2.append(h)
                S_3.append(h)
                S_4.append(h - (y_inter[i + 1] - y_inter[i]))

            if (x_inter[i] != x_integ[i] and x_inter[j] == x_integ[j] and
                        y_inter[i] == y_integ[i] and y_inter[j] != y_integ[j]):
                #                   CASE 15
                case.append(15)
                x_case15.append(x_integ[i] + h / 2)
                y_case15.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] + h / 2)
                V.append(0.5 * (x_inter[i + 1] - x_inter[i]) *
                         (y_inter[i + 1] - y_inter[i]))
                S.append(np.sqrt((x_inter[i + 1] - x_inter[i]) ** 2 +
                                 (y_inter[i + 1] - y_inter[i]) ** 2))

                S_1.append(0)
                S_2.append(y_inter[i + 1] - y_inter[i])
                S_3.append(x_inter[i + 1] - x_inter[i])
                S_4.append(0)

            if (x_inter[i] != x_integ[i] and x_inter[j] != x_integ[j] and
                        y_inter[i] == y_integ[i] and y_inter[j] == y_integ[j]):
                #                   CASE 16
                case.append(16)
                x_case16.append(x_integ[i] + h / 2)
                y_case16.append(y_integ[i] + h / 2)
                x_cc.append(x_integ[i] + h / 2)
                y_cc.append(y_integ[i] + h / 2)
                V.append(h ** 2 - (0.5 * h * (x_inter[i + 1] - x_inter[i]) +
                                   h * (x_inter[i] - x_integ[i])))
                S.append(np.sqrt(h ** 2 + (x_inter[i + 1] - x_inter[i]) ** 2))

                S_1.append(h - (x_inter[i + 1] - x_integ[i + 1]))
                S_2.append(h)
                S_3.append(h - (x_inter[i] - x_integ[i]))
                S_4.append(0)

    x_halcell = x_case1 + x_case2 + x_case4 + \
                x_case5 + x_case6 + x_case8 + \
                x_case9 + x_case10 + x_case12 + \
                x_case13 + x_case14 + x_case16

    y_halcell = y_case1 + y_case2 + y_case4 + \
                y_case5 + y_case6 + y_case8 + \
                y_case9 + y_case10 + y_case12 + \
                y_case13 + y_case14 + y_case16

    x_cutcell = x_case1 + x_case2 + x_case3 + x_case4 + \
                x_case5 + x_case6 + x_case7 + x_case8 + \
                x_case9 + x_case10 + x_case11 + x_case12 + \
                x_case13 + x_case14 + x_case15 + x_case16

    y_cutcell = y_case1 + y_case2 + y_case3 + y_case4 + \
                y_case5 + y_case6 + y_case7 + y_case8 + \
                y_case9 + y_case10 + y_case11 + y_case12 + \
                y_case13 + y_case14 + y_case15 + y_case16

    plot_coordinates=False
    if plot_coordinates:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
        fig = plt.figure('GRID vs. CIRCLE (2-B)')
        ax1 = fig.add_subplot(111)
        ax1.set_aspect('equal')

        plt.plot(ax, by, 'ko', markersize=10, markeredgewidth=4, markeredgecolor='k', markerfacecolor='None')
        #    plt.plot(x_circ,y_circ,'c-')
        # plt.plot(x_cutcell,y_cutcell,'ro', markersize=3, markeredgewidth=4, markeredgecolor='r')
        # plt.plot(x_halcell, y_halcell, 'gs', markersize=3, markeredgewidth=4, markeredgecolor='g')

        plt.plot(x_out, y_out, 'md--')

        ax1.xaxis.set_ticks(np.arange(0, 2 * ax + h, h))
        ax1.yaxis.set_ticks(np.arange(0, 2 * by + h, h))
        plt.grid(True, which='major', color='b', linestyle='--', lw=1.)

        plt.minorticks_on()
        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=3.)
        plt.xlim(ax - 1.2 * radius, ax + 1.4 * radius)
        plt.ylim(by - 1.2 * radius, by + 1.4 * radius)
        plt.show()

    return x_inter, y_inter, \
           x_integ, y_integ, \
           x_ins, y_ins, \
           x_out, y_out, \
           case, \
           h, radius, \
           x_halcell, y_halcell, \
           x_cutcell, y_cutcell, \
           x_cc, y_cc, \
           S_1, S_2, S_3, S_4, \
           S, V