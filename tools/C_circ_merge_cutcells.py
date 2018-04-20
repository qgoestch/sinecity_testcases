# -*- coding: utf-8 -*-
##
# \file     C_circ_merge_cutcells.py
# \title    Get the parameters of the cut cells after merging them.
#           This parameters are required for FVTD calculations only, but is mandatory to get the corner
#           of the FDTD circle.
# \author   Pierre Chobeau
# \version  0.1
# \license  BSD 3-Clause License
# \inst     UMRAE (Ifsttar Nantes), LAUM (Le Mans Université)
# \date     2017, 07 Sep.
##

import numpy as np
import B_circ_cut_cell_types as B_5

def merging_cutcells(radius, h, ax, by):
    """
          Step 3/4, calculation of the cut-cells parameters - after merging - for the FVTD update.

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

    :return     the new parameters (volumes, surfaces, corrdinates) of the cut cells after merging.
    """
    x_inter, y_inter, \
    x_integ, y_integ, \
    x_ins, y_ins, \
    x_out, y_out, \
    case, \
    h, radius, \
    x_halcell, y_halcell, \
    x_cutcell, y_cutcell, \
    x_cc, y_cc, \
    S_1, S_2, S_3, S_4, \
    S, V = B_5.cut_cell_types(radius, h, ax, by)

    # ==============================================================================
    #                       MERGING OF THE CELLS
    #   - CALCULATION OF THE MERGED VOLUME AND SURFACES FOR EACH MERGED CELL;
    #   - CALCULATION OF THE DISTANCES h_1, h_2, h_3, h_4; <-- !!! NO NEED !!!
    #   - CALCULATION OF THE CENTERS OF THE MERGED CELLS x_cc, y_cc.
    # ==============================================================================
    area_x_r = [];
    area_y_r = []
    area_x_y = [];
    area_y_y = []
    area_x_g = [];
    area_y_g = []
    area_x_m = [];
    area_y_m = []
    area_x_c = [];
    area_y_c = []
    area_x_b = [];
    area_y_b = []
    area_x_k = [];
    area_y_k = []
    area_x_gray = [];
    area_y_gray = []

    V_b_m = [];
    S_b_m = [];
    x_cc_m = [];
    y_cc_m = []
    S_1_m = [];
    S_2_m = [];
    S_3_m = [];
    S_4_m = [];
    S_5_m = [];
    S_6_m = []
    case_color = []
    x_cc_tri = [];
    y_cc_tri = []
    S_1_tri = [];
    S_2_tri = [];
    S_3_tri = [];
    S_4_tri = []

    blue_no_pb = True
    for i in range(len(V)):
        #       1st quadrant
        if ((y_integ[i] - by) > 0 and (x_integ[i] - ax) > 0):
            if i == len(x_inter) - 1:
                j = 0
            else:
                j = i + 1

            if case[i] == 1 and case[i + 1] == 1:
                # ==============================================================================
                #               CASE 1 = RED
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i] + h)
                    S_2_m.append(S_2[i])
                    S_3_m.append(S_3[i] + h)
                    S_4_m.append(S_4[i])
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])

                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i] + h)
                    y_cc_m.append(y_cc[i])

                    area_x_r.append([x_inter[i], x_inter[i + 1], x_out[i + 1] + h, x_out[i] + h])
                    area_y_r.append([y_inter[i], y_inter[i + 1], y_out[i + 1], y_out[i]])
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(1)

            if case[i] == 1 and case[i - 1] == 1 and case[i + 1] == 2:
                # ==============================================================================
                #           CASE 2 = YELLOW
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(h)
                S_2_m.append(h)
                S_3_m.append(S_3[i] + h)
                S_4_m.append(0)
                S_5_m.append(0)
                S_6_m.append(S6)

                S_b_m.append(S[i] - Sb_spt6)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt6 * S6 / 2))

                x_cc_m.append(x_cc[i] + h)
                y_cc_m.append(y_cc[i])

                area_x_y.append([x_inter[i], xc6, x_out[i + 1],
                                 x_out[i + 1] + h, x_out[i] + h])
                area_y_y.append([y_inter[i], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i]])
                area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])

                case_color.append(2)

            if (case[i] == 2 and (case[i - 1] == 1 and case[i - 2] != 1)
                and (case[i + 1] == 3 or case[i + 1] == 4)):
                # ==============================================================================
                #           CASE 3 = BLUE1
                # ==============================================================================
                blue_no_pb = False
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                A5 = x_out[i - 1] - ax
                C5 = np.sqrt((x_out[i - 1] - ax) ** 2 + (y_out[i - 1] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i - 1] - xc5) ** 2 +
                             (y_out[i - 1] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i - 1] - xc5) ** 2 +
                                  (y_inter[i - 1] - yc5) ** 2)
                S_1_m.append(h)
                S_2_m.append(2. * h)
                S_3_m.append(0)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i - 1] + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i - 1] + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1] + h, x_out[i - 1]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i - 1]])
                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(3)

            if (
                            (case[i] == 2 and (case[i - 2] == 2 and case[i - 1] == 3 and
                                                       case[i + 1] == 3 and case[i + 2] == 2)) or
                            (case[i] == 2 and (case[i - 2] == 1 and case[i - 1] == 1 and
                                                       case[i + 1] == 3 and (case[i + 2] == 1 or case[i + 2] == 2))) or
                        (case[i] == 2 and ((case[i - 2] == 4 or case[i - 2] == 2) and case[i - 1] == 3 and
                                                   case[i + 1] == 4 and case[i + 2] == 4))
            ):
                # ==============================================================================
                #           CASE 4 = MAGENTA
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(h)
                S_2_m.append(h)
                S_3_m.append(0)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + (Sb_spt6 * S6 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])

                area_x_m.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1] + h, x_out[i]])
                area_y_m.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i]])

                if case[i - 1] == 1 or case[i + 1] == 4:
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i] + h, y_out[i]])
                else:
                    area_x_gray.append([x_out[i], x_out[i + 1], x_out[i]])
                    area_y_gray.append([y_out[i] + h, y_out[i + 1], y_out[i]])
                    area_x_b.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i], x_out[i]])
                    area_y_b.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i] + h, y_out[i]])

                case_color.append(4)

            if (case[i] == 2 and (case[i + 1] == 4 and case[i + 2] != 4)
                and (case[i + 1] == 3 or case[i + 1] == 4)):
                # ==============================================================================
                #           CASE 5 = BLUE2
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 2] - ax
                C6 = np.sqrt((x_out[i + 2] - ax) ** 2 + (y_out[i + 2] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 2] - xc6) ** 2 +
                             (y_out[i + 2] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 2] - xc6) ** 2 +
                                  (y_inter[i + 2] - yc6) ** 2)

                S_1_m.append(2. * h)
                S_2_m.append(h)
                S_3_m.append(0)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + S[i + 1] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + V[i + 1] + (Sb_spt6 * S6 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 2],
                                 x_out[i + 1] + h, x_out[i]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 2],
                                 y_out[i + 1], y_out[i]])
                case_color.append(5)

            if case[i] == 4 and case[i + 1] == 4 and case[i - 1] != 4:
                # ==============================================================================
                #           CASE 6 = GREEN
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                S_1_m.append(h)
                S_2_m.append(h)
                S_3_m.append(0)
                S_4_m.append(S_4[i] + h)
                S_5_m.append(S5)
                S_6_m.append(0)
                S_b_m.append(S[i] - Sb_spt5)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt5 * S5 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i] + h)

                area_x_g.append([xc5, x_inter[i + 1], x_inter[i + 1],
                                 x_out[i + 1], x_out[i], x_out[i]])
                area_y_g.append([yc5, y_inter[i + 1], y_inter[i + 1],
                                 y_out[i + 1] + h, y_out[i] + h, y_out[i]])
                area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])

                case_color.append(6)

            if case[i] == 4 and case[i - 1] == 4:
                # ==============================================================================
                #               CASE 7 = CYAN
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i])
                    S_2_m.append(S_2[i] + h)
                    S_3_m.append(S_3[i])
                    S_4_m.append(S_4[i] + h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i])
                    y_cc_m.append(y_cc[i] + h)

                    area_x_c.append([x_inter[i], x_inter[i + 1], x_out[i + 1], x_out[i]])
                    area_y_c.append([y_inter[i], y_inter[i + 1], y_out[i + 1] + h, y_out[i] + h])
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(7)

                    #       2nd quadrant
        if ((y_integ[i] - by) > 0 and (x_integ[i] - ax) < 0):
            if case[i] == 5 and case[i + 1] == 5:
                # ==============================================================================
                #               CASE 8 = RED
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i])
                    S_2_m.append(S_2[i] + h)
                    S_3_m.append(S_3[i])
                    S_4_m.append(S_4[i] + h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i])
                    y_cc_m.append(y_cc[i] + h)

                    area_x_r.append([x_inter[i], x_inter[i + 1], x_out[i + 1], x_out[i]])
                    area_y_r.append([y_inter[i], y_inter[i + 1], y_out[i + 1] + h, y_out[i] + h])
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(8)

            if case[i] == 5 and case[i - 1] == 5 and case[i + 1] == 6:
                # ==============================================================================
                #           CASE 9 = YELLOW
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(h)
                S_2_m.append(S_2[i] + h)
                S_3_m.append(0)
                S_4_m.append(h)
                S_5_m.append(0)
                S_6_m.append(S6)

                S_b_m.append(S[i] - Sb_spt6)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt6 * S6 / 2))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i] + h)

                area_x_y.append([x_inter[i], xc6, x_out[i + 1],
                                 x_out[i + 1], x_out[i]])
                area_y_y.append([y_inter[i], yc6, y_out[i + 1],
                                 y_out[i + 1] + h, y_out[i] + h])
                case_color.append(9)

            if (case[i] == 6 and (case[i - 1] == 5 and case[i - 2] != 5)
                and (case[i + 1] == 7 or case[i + 1] == 8)):
                # ==============================================================================
                #           CASE 10 = BLUE1
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                A5 = x_out[i - 1] - ax
                C5 = np.sqrt((x_out[i - 1] - ax) ** 2 + (y_out[i - 1] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i - 1] - xc5) ** 2 +
                             (y_out[i - 1] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i - 1] - xc5) ** 2 +
                                  (y_inter[i - 1] - yc5) ** 2)
                S_1_m.append(2. * h)
                S_2_m.append(0)
                S_3_m.append(0)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i - 1] + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i - 1] + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1], x_out[i - 1]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1] + h, y_out[i - 1]])
                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(10)

            if (
                            (case[i] == 6 and (case[i - 2] == 6 and case[i - 1] == 7 and
                                                       case[i + 1] == 7 and case[i + 2] == 6)) or
                            (case[i] == 6 and (case[i - 2] == 5 and case[i - 1] == 5 and
                                                       case[i + 1] == 7 and (case[i + 2] == 5 or case[i + 2] == 6))) or
                        (case[i] == 6 and ((case[i - 2] == 8 or case[i - 2] == 6) and case[i - 1] == 7 and
                                                   case[i + 1] == 8 and case[i + 2] == 8))
            ):
                # ==============================================================================
                #           CASE 11 = MAGENTA
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(h)
                S_2_m.append(0)
                S_3_m.append(0)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_m.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i] - h, x_out[i]])
                area_y_m.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i], y_out[i]])

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(11)

            if (case[i] == 6 and (case[i + 1] == 8 and case[i + 2] != 8)
                and (case[i + 1] == 7 or case[i + 1] == 8)):
                # ==============================================================================
                #           CASE 12 = BLUE2
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 2] - ax
                C6 = np.sqrt((x_out[i + 2] - ax) ** 2 + (y_out[i + 2] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 2] - xc6) ** 2 +
                             (y_out[i + 2] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 2] - xc6) ** 2 +
                                  (y_inter[i + 2] - yc6) ** 2)

                S_1_m.append(h)
                S_2_m.append(0)
                S_3_m.append(0)
                S_4_m.append(2. * h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + S[i + 1] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + V[i + 1] + (Sb_spt6 * S6 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 2],
                                 x_out[i] - h, x_out[i]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 2],
                                 y_out[i], y_out[i]])
                case_color.append(12)

            if case[i] == 8 and case[i + 1] == 8 and case[i - 1] != 8:
                # ==============================================================================
                #           CASE 13 = GREEN
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                S_1_m.append(h)
                S_2_m.append(0)
                S_3_m.append(S_3[i] + h)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(0)
                S_b_m.append(S[i] - Sb_spt5)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt5 * S5 / 2.))

                x_cc_m.append(x_cc[i] - h)
                y_cc_m.append(y_cc[i])

                area_x_g.append([xc5, x_inter[i + 1],
                                 x_out[i + 1] - h, x_out[i] - h, x_out[i]])
                area_y_g.append([yc5, y_inter[i + 1],
                                 y_out[i + 1], y_out[i], y_out[i]])
                case_color.append(13)

            if case[i] == 8 and case[i - 1] == 8:
                # ==============================================================================
                #               CASE 14 = CYAN
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i] + h)
                    S_2_m.append(0)
                    S_3_m.append(S_3[i] + h)
                    S_4_m.append(h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i] - h)
                    y_cc_m.append(y_cc[i])

                    area_x_c.append([x_inter[i], x_inter[i + 1], x_out[i + 1] - h, x_out[i] - h])
                    area_y_c.append([y_inter[i], y_inter[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(14)

                    #       3rd quadrant
        if ((y_integ[i] - by) < 0 and (x_integ[i] - ax) < 0):
            if case[i] == 9 and case[i + 1] == 9:
                # ==============================================================================
                #               CASE 15 = RED
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i] + h)
                    S_2_m.append(0)
                    S_3_m.append(S_3[i] + h)
                    S_4_m.append(h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i] - h)
                    y_cc_m.append(y_cc[i])

                    area_x_r.append([x_inter[i], x_inter[i + 1], x_out[i + 1] - h, x_out[i] - h])
                    area_y_r.append([y_inter[i], y_inter[i + 1], y_out[i + 1], y_out[i]])
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(15)

            if case[i] == 9 and case[i - 1] == 9 and case[i + 1] == 10:
                # ==============================================================================
                #           CASE 16 = YELLOW
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = - radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(S_1[i] + h)
                S_2_m.append(0)
                S_3_m.append(h)
                S_4_m.append(h)
                S_5_m.append(0)
                S_6_m.append(S6)

                S_b_m.append(S[i] - Sb_spt6)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt6 * S6 / 2))

                x_cc_m.append(x_cc[i] - h)
                y_cc_m.append(y_cc[i])

                area_x_y.append([x_inter[i], xc6, x_out[i + 1],
                                 x_out[i + 1] - h, x_out[i] - h])
                area_y_y.append([y_inter[i], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i]])
                case_color.append(16)

            if (case[i] == 10 and (case[i - 1] == 9 and case[i - 2] != 9)
                and (case[i + 1] == 11 or case[i + 1] == 12)):
                # ==============================================================================
                #           CASE 17 = BLUE1
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                A5 = x_out[i - 1] - ax
                C5 = np.sqrt((x_out[i - 1] - ax) ** 2 + (y_out[i - 1] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i - 1] - xc5) ** 2 +
                             (y_out[i - 1] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i - 1] - xc5) ** 2 +
                                  (y_inter[i - 1] - yc5) ** 2)
                S_1_m.append(0)
                S_2_m.append(0)
                S_3_m.append(h)
                S_4_m.append(2. * h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i - 1] + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i - 1] + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1] - h, x_out[i - 1]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i - 1]])
                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(17)

            if (
                            (case[i] == 10 and (case[i - 2] == 10 and case[i - 1] == 11 and
                                                        case[i + 1] == 11 and case[i + 2] == 10)) or
                            (case[i] == 10 and (case[i - 2] == 9 and case[i - 1] == 9 and
                                                        case[i + 1] == 11 and (
                                    case[i + 2] == 9 or case[i + 2] == 10))) or
                        (case[i] == 10 and ((case[i - 2] == 12 or case[i - 2] == 10) and case[i - 1] == 11 and
                                                    case[i + 1] == 12 and case[i + 2] == 12))
            ):
                # ==============================================================================
                #           CASE 18 = MAGENTA
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(0)
                S_2_m.append(0)
                S_3_m.append(h)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_m.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1] - h, x_out[i]])
                area_y_m.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1], y_out[i]])

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(18)

            if (case[i] == 10 and (case[i + 1] == 12 and case[i + 2] != 12)
                and (case[i + 1] == 11 or case[i + 1] == 12)):
                # ==============================================================================
                #           CASE 19 = BLUE2
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 2] - ax
                C6 = np.sqrt((x_out[i + 2] - ax) ** 2 + (y_out[i + 2] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 2] - xc6) ** 2 +
                             (y_out[i + 2] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 2] - xc6) ** 2 +
                                  (y_inter[i + 2] - yc6) ** 2)

                S_1_m.append(0)
                S_2_m.append(0)
                S_3_m.append(2. * h)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + S[i + 1] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + V[i + 1] + (Sb_spt6 * S6 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 2],
                                 x_out[i], x_out[i]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 2],
                                 y_out[i] - h, y_out[i]])
                case_color.append(19)

            if case[i] == 12 and case[i + 1] == 12 and case[i - 1] != 12:
                # ==============================================================================
                #           CASE 20 = GREEN
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                S_1_m.append(0)
                S_2_m.append(S_2[i] + h)
                S_3_m.append(h)
                S_4_m.append(h)
                S_5_m.append(S5)
                S_6_m.append(0)
                S_b_m.append(S[i] - Sb_spt5)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt5 * S5 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i] - h)

                area_x_g.append([xc5, x_inter[i + 1],
                                 x_out[i + 1], x_out[i], x_out[i]])
                area_y_g.append([yc5, y_inter[i + 1],
                                 y_out[i + 1] - h, y_out[i] - h, y_out[i]])
                case_color.append(20)

            if case[i] == 12 and case[i - 1] == 12:
                # ==============================================================================
                #               CASE 21 = CYAN
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(0)
                    S_2_m.append(S_2[i] + h)
                    S_3_m.append(h)
                    S_4_m.append(S_4[i] + h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i])
                    y_cc_m.append(y_cc[i] - h)

                    area_x_c.append([x_inter[i], x_inter[i + 1], x_out[i + 1], x_out[i]])
                    area_y_c.append([y_inter[i], y_inter[i + 1], y_out[i + 1] - h, y_out[i] - h])
                    case_color.append(21)


                    #       4th quadrant
        if ((y_integ[i] - by) < 0 and (x_integ[i] - ax) > 0):
            j = i + 1;
            k = i + 2
            if i == len(case) - 1:    j = 0; k = 1;
            if i == len(case) - 2:    k = 0
            if case[i] == 13 and case[i + 1] == 13:
                # ==============================================================================
                #               CASE 22 = RED
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(0)
                    S_2_m.append(S_2[i] + h)
                    S_3_m.append(h)
                    S_4_m.append(S_4[i] + h)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i])
                    y_cc_m.append(y_cc[i] - h)

                    area_x_r.append([x_inter[i], x_inter[i + 1], x_out[i + 1], x_out[i]])
                    area_y_r.append([y_inter[i], y_inter[i + 1], y_out[i + 1] - h, y_out[i] - h])
                    case_color.append(22)

            if case[i] == 13 and case[i - 1] == 13 and case[i + 1] == 14:
                # ==============================================================================
                #           CASE 23 = YELLOW
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = - radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                S_1_m.append(0)
                S_2_m.append(h)
                S_3_m.append(h)
                S_4_m.append(S_4[i] + h)
                S_5_m.append(0)
                S_6_m.append(S6)

                S_b_m.append(S[i] - Sb_spt6)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt6 * S6 / 2))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i] - h)

                area_x_y.append([x_inter[i], xc6, x_out[i + 1],
                                 x_out[i + 1], x_out[i]])
                area_y_y.append([y_inter[i], yc6, y_out[i + 1],
                                 y_out[i + 1] - h, y_out[i] - h])
                case_color.append(23)

            if (case[i] == 14 and (case[i - 1] == 13 and case[i - 2] != 13)
                and (case[i + 1] == 15 or case[i + 1] == 16)):
                # ==============================================================================
                #           CASE 24 = BLUE1
                # ==============================================================================
                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[i + 1] - ax) ** 2 + (y_out[i + 1] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[i + 1] - xc6) ** 2 +
                             (y_out[i + 1] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[i + 1] - xc6) ** 2 +
                                  (y_inter[i + 1] - yc6) ** 2)

                A5 = x_out[i - 1] - ax
                C5 = np.sqrt((x_out[i - 1] - ax) ** 2 + (y_out[i - 1] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i - 1] - xc5) ** 2 +
                             (y_out[i - 1] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i - 1] - xc5) ** 2 +
                                  (y_inter[i - 1] - yc5) ** 2)
                S_1_m.append(0)
                S_2_m.append(h)
                S_3_m.append(2. * h)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i - 1] + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i - 1] + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_b.append([xc5, x_inter[i], x_inter[i + 1], xc6, x_out[i + 1],
                                 x_out[i + 1], x_out[i - 1]])
                area_y_b.append([yc5, y_inter[i], y_inter[i + 1], yc6, y_out[i + 1],
                                 y_out[i + 1] - h, y_out[i - 1]])
                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(24)

            if (
                            (case[i] == 14 and (case[i - 2] == 14 and case[i - 1] == 15 and
                                                        case[j] == 15 and case[k] == 14)) or
                            (case[i] == 14 and (case[i - 2] == 13 and case[i - 1] == 13 and
                                                        case[j] == 15 and (case[k] == 13 or case[k] == 14))) or
                        (case[i] == 14 and ((case[i - 2] == 16 or case[i - 2] == 14) and case[i - 1] == 15 and
                                                    case[j] == 16 and case[k] == 16))
            ):
                # ==============================================================================
                #           CASE 25 = MAGENTA
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[i + 1] - ax
                C6 = np.sqrt((x_out[j] - ax) ** 2 + (y_out[j] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[j] - xc6) ** 2 +
                             (y_out[j] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[j] - xc6) ** 2 +
                                  (y_inter[j] - yc6) ** 2)

                S_1_m.append(0)
                S_2_m.append(h)
                S_3_m.append(h)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + (Sb_spt6 * S6 / 2.))

                area_x_m.append([xc5, x_inter[i], x_inter[j], xc6, x_out[j],
                                 x_out[i] + h, x_out[i]])
                area_y_m.append([yc5, y_inter[i], y_inter[j], yc6, y_out[j],
                                 y_out[i], y_out[i]])

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])
                case_color.append(25)

            if (case[i] == 14 and (case[j] == 16 and case[k] != 16)
                and (case[j] == 15 or case[j] == 16)):
                # ==============================================================================
                #           CASE 26 = BLUE2
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                A6 = x_out[k] - ax
                C6 = np.sqrt((x_out[k] - ax) ** 2 + (y_out[k] - by) ** 2)
                Al6 = np.arccos(A6 / C6)
                xc6 = radius * np.cos(Al6) + ax
                yc6 = -radius * np.sin(Al6) + by
                S6 = np.sqrt((x_out[k] - xc6) ** 2 +
                             (y_out[k] - yc6) ** 2)
                Sb_spt6 = np.sqrt((x_inter[k] - xc6) ** 2 +
                                  (y_inter[k] - yc6) ** 2)

                S_1_m.append(0)
                S_2_m.append(2. * h)
                S_3_m.append(h)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(S6)
                S_b_m.append(Sb_spt5 + S[i] + S[j] + Sb_spt6)
                V_b_m.append((Sb_spt5 * S5 / 2.) + V[i] + V[j] + (Sb_spt6 * S6 / 2.))

                x_cc_m.append(x_cc[i])
                y_cc_m.append(y_cc[i])

                area_x_b.append([xc5, x_inter[i], x_inter[j], xc6, x_out[k],
                                 x_out[i] + h, x_out[i]])
                area_y_b.append([yc5, y_inter[i], y_inter[j], yc6, y_out[k],
                                 y_out[i], y_out[i]])
                case_color.append(26)

            if case[i] == 16 and case[j] == 16 and case[i - 1] != 16:
                # ==============================================================================
                #           CASE 27 = GREEN
                # ==============================================================================
                A5 = x_out[i + 0] - ax
                C5 = np.sqrt((x_out[i + 0] - ax) ** 2 + (y_out[i + 0] - by) ** 2)
                Al5 = np.arccos(A5 / C5)
                xc5 = radius * np.cos(Al5) + ax
                yc5 = -radius * np.sin(Al5) + by
                S5 = np.sqrt((x_out[i + 0] - xc5) ** 2 +
                             (y_out[i + 0] - yc5) ** 2)
                Sb_spt5 = np.sqrt((x_inter[i + 0] - xc5) ** 2 +
                                  (y_inter[i + 0] - yc5) ** 2)

                S_1_m.append(S_1[i] + h)
                S_2_m.append(h)
                S_3_m.append(h)
                S_4_m.append(0)
                S_5_m.append(S5)
                S_6_m.append(0)
                S_b_m.append(S[i] - Sb_spt5)
                V_b_m.append(V[i] + h ** 2 - (Sb_spt5 * S5 / 2.))

                x_cc_m.append(x_cc[i] + h)
                y_cc_m.append(y_cc[i])

                area_x_g.append([xc5, x_inter[j],
                                 x_out[j] + h, x_out[i] + h, x_out[i]])
                area_y_g.append([yc5, y_inter[j],
                                 y_out[j], y_out[i], y_out[i]])
                case_color.append(27)

            if case[i] == 16 and case[i - 1] == 16:
                # ==============================================================================
                #               CASE 28 = CYAN
                # ==============================================================================
                if V[i] < (.9 * h ** 2):
                    S_1_m.append(S_1[i] + h)
                    S_2_m.append(h)
                    S_3_m.append(S_3[i] + h)
                    S_4_m.append(0)
                    S_5_m.append(0)
                    S_6_m.append(0)

                    S_b_m.append(S[i])
                    V_b_m.append(V[i] + h ** 2)

                    x_cc_m.append(x_cc[i] + h)
                    y_cc_m.append(y_cc[i])

                    area_x_c.append([x_inter[i], x_inter[i + 1], x_out[i + 1] + h, x_out[i] + h])
                    area_y_c.append([y_inter[i], y_inter[i + 1], y_out[i + 1], y_out[i]])
                    area_x_k.append([x_ins[i], x_ins[i + 1], x_out[i + 1], x_out[i]])
                    area_y_k.append([y_ins[i], y_ins[i + 1], y_out[i + 1], y_out[i]])
                    case_color.append(28)

    if blue_no_pb == True:
        print h

    # ==============================================================================
    #   Plot the 1st quadrant for staircases, halved corner and cut cells.
    #   Figure 1 of the rejected JASAEL --> tech. rep. of Nov. 2016
    # ==============================================================================
    plot = False
    if plot:
        phi = np.linspace(0, 2 * np.pi, 10000)
        x_circ = np.zeros((len(phi)), dtype=np.float64)
        x_circ = radius * np.cos(phi) + ax
        y_circ = np.zeros((len(phi)), dtype=np.float64)
        y_circ = radius * np.sin(phi) + by

        import matplotlib.pyplot as plt
        import pylab as pyl
        fig = plt.figure('Meshes', figsize=(10, 4.1))
        #        fig = plt.figure('sc')
        ax1 = fig.add_subplot(131)
        ax1.set_aspect('equal')

        ax1.xaxis.set_ticks(np.arange(0, 2 * ax + h, h))
        ax1.yaxis.set_ticks(np.arange(0, 2 * by + h, h))

        plt.minorticks_on()
        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=.3)

        for i in range(len(area_x_k)):
            pyl.fill(area_x_k[i], area_y_k[i], 'k', alpha=0.8, edgecolor='k')

        for i in range(len(area_x_b)):
            pyl.fill(area_x_b[i], area_y_b[i], 'k', alpha=0.8, edgecolor='k')

        plt.plot(ax, by, 'kx', markersize=5, markeredgewidth=2, markeredgecolor='k',
                 markerfacecolor='None')
        plt.plot(x_circ, y_circ, 'r-', lw=2)
        plt.xticks(np.arange(ax, ax + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.yticks(np.arange(by, by + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.xlabel("Cell number", fontsize=14)
        plt.ylabel("Cell number", fontsize=14)
        plt.xlim(ax - 1.2 * h, ax + 1.2 * radius)
        plt.ylim(by - 1.2 * h, by + 1.2 * radius)

        #        fig = plt.figure('hc')
        ax1 = fig.add_subplot(132)
        ax1.set_aspect('equal')

        ax1.xaxis.set_ticks(np.arange(0, 2 * ax + h, h))
        ax1.yaxis.set_ticks(np.arange(0, 2 * by + h, h))

        plt.minorticks_on()
        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=.3)

        for i in range(len(area_x_k)):
            pyl.fill(area_x_k[i], area_y_k[i], 'k', alpha=0.8, edgecolor='k')

        for i in range(len(area_x_gray)):
            pyl.fill(area_x_gray[i], area_y_gray[i], facecolor='0.65', edgecolor='k', linewidth='1.5')

        plt.plot(ax, by, 'kx', markersize=5, markeredgewidth=2, markeredgecolor='k',
                 markerfacecolor='None')
        plt.plot(x_circ, y_circ, 'r-', lw=2)
        plt.xticks(np.arange(ax, ax + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.yticks(np.arange(by, by + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.xlabel("Cell number", fontsize=14)
        plt.xlim(ax - 1.2 * h, ax + 1.2 * radius)
        plt.ylim(by - 1.2 * h, by + 1.2 * radius)

        #        fig = plt.figure('cc')
        ax1 = fig.add_subplot(133)
        ax1.set_aspect('equal')

        #    plt.plot(xc,yc,'bo')
        #        plt.plot(x_inter,y_inter,'ro')
        #        plt.plot(x_ins,y_ins,'cd')
        #        plt.plot(x_out,y_out,'ys')
        #    plt.plot(x_integ,y_integ,'bs')

        ax1.xaxis.set_ticks(np.arange(0, 2 * ax + h, h))
        ax1.yaxis.set_ticks(np.arange(0, 2 * by + h, h))

        plt.minorticks_on()
        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=.3)

        for i in range(len(area_x_r)):
            pyl.fill(area_x_r[i], area_y_r[i], facecolor='0.65', edgecolor='k', linewidth='1.5')

        for i in range(len(area_x_y)):
            pyl.fill(area_x_y[i], area_y_y[i], facecolor='0.65', edgecolor='k', linewidth='1.5')

        for i in range(len(area_x_g)):
            pyl.fill(area_x_g[i], area_y_g[i], facecolor='0.65', edgecolor='k', linewidth='1.5')

        for i in range(len(area_x_m)):
            pyl.fill(area_x_m[i], area_y_m[i], facecolor='0.65', edgecolor='k', linewidth='1.5')

        for i in range(len(area_x_c)):
            pyl.fill(area_x_c[i], area_y_c[i], facecolor='0.65', edgecolor='k', linewidth='1.5')
        # for i in range(len(area_x_r)):
        #            pyl.fill(area_x_r[i],area_y_r[i], 'r', alpha=0.7, edgecolor='k')
        #
        #        for i in range(len(area_x_y)):
        #            pyl.fill(area_x_y[i],area_y_y[i], 'y', alpha=0.7, edgecolor='k')
        #
        #        for i in range(len(area_x_g)):
        #            pyl.fill(area_x_g[i],area_y_g[i], 'g', alpha=0.7, edgecolor='k')
        #
        #        for i in range(len(area_x_m)):
        #            pyl.fill(area_x_m[i],area_y_m[i], 'm', alpha=0.7, edgecolor='k')
        #
        #        for i in range(len(area_x_c)):
        #            pyl.fill(area_x_c[i],area_y_c[i], 'c', alpha=0.7, edgecolor='k')

        #        for i in range(len(area_x_b)):
        #            pyl.fill(area_x_b[i],area_y_b[i], 'b', alpha=0.7, edgecolor='k')

        #        plt.plot(x_cc_m[:], y_cc_m[:], 'ko', markersize=3, markeredgewidth=2,
        #                 markeredgecolor='k', markerfacecolor='None')

        plt.plot(ax, by, 'kx', markersize=5, markeredgewidth=2, markeredgecolor='k',
                 markerfacecolor='None')
        plt.plot(x_circ, y_circ, 'r-', lw=2)
        plt.xticks(np.arange(ax, ax + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.yticks(np.arange(by, by + 12 * h, h), ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
        plt.xlabel("Cell number", fontsize=14)
        plt.xlim(ax - 1.2 * h, ax + 1.2 * radius)
        plt.ylim(by - 1.2 * h, by + 1.2 * radius)

    # plt.show()
    #        plt.savefig('Figure1.jpeg',transparent=True,bbox_inches='tight',pad_inches = 0)
    #        plt.savefig('Figure1.eps',transparent=True,bbox_inches='tight',pad_inches = 0)
    #        plt.savefig('Figure1.pdf',transparent=True,bbox_inches='tight',pad_inches = 0)

    # ==============================================================================
    #   ROUNDING FOR MATCHING THE SURFACES BEFORE THE TEST
    # ==============================================================================
    S_1_m = [round(i, 4) for i in S_1_m]
    S_2_m = [round(i, 4) for i in S_2_m]
    S_3_m = [round(i, 4) for i in S_3_m]
    S_4_m = [round(i, 4) for i in S_4_m]
    S_1_tri = [round(i, 4) for i in S_1_tri]
    S_2_tri = [round(i, 4) for i in S_2_tri]
    S_3_tri = [round(i, 4) for i in S_3_tri]
    S_4_tri = [round(i, 4) for i in S_4_tri]
    # ==============================================================================
    #           TEST OF THE CONECTED SURFACES BTW BOUNDARY-CELLS
    #               DOES THE BOUNDARY-CELLS SURFACES MATCH?
    #         TEST OK FOR h_step=[0.0355,0.0405,0.0475,0.0550,0.0675]
    # ==============================================================================
    test = False
    if test:
        for i in range(len(case_color)):
            if case_color[i] == 1:
                print 'case_c=%i' % (case_color[i])
                if S_2_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_1_m[i] != S_3_m[i + 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i + 1])), 20)
                if S_3_m[i] != S_1_m[i - 1]: print 'error on S_3_m: %.20f' % round((abs(S_3_m[i] - S_1_m[i - 1])), 20)
            if case_color[i] == 2:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m=%f and h=%f' % (S_1_m[i], h)
                if S_2_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_3_m[i] != S_1_m[i - 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i + 1])), 20)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
            if case_color[i] == 3:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m'
                if S_2_m[i] != 2. * h: print 'error on S_2_m=%.20f and h=%.20f' % (S_2_m[i], h)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 4:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m'
                if S_2_m[i] != h: print 'error on S_2_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 5:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != 2. * h: print 'error on S_1_m'
                if S_2_m[i] != h: print 'error on S_2_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 6:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m=%f and h=%f' % (S_1_m[i], h)
                if S_2_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_4_m[i] != S_2_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i + 1])), 4)
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 7:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m'
                if S_4_m[i] != S_2_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i + 1])), 4)
                if S_2_m[i] != S_4_m[i - 1]: print 'error on S_2_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i - 1])), 4)

            if case_color[i] == 8:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m'
                if S_4_m[i] != S_2_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i + 1])), 20)
                if S_2_m[i] != S_4_m[i - 1]: print 'error on S_2_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i - 1])), 20)
            if case_color[i] == 9:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m=%f and h=%f' % (S_1_m[i], h)
                if S_4_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_2_m[i] != S_4_m[i - 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i + 1])), 20)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
            if case_color[i] == 10:
                print 'case_c=%i' % (case_color[i])
                if S_4_m[i] != h: print 'error on S_1_m'
                if S_1_m[i] != 2. * h: print 'error on S_2_m=%.20f and h=%.20f' % (S_2_m[i], h)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 11:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m'
                if S_4_m[i] != h: print 'error on S_2_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 12:
                print 'case_c=%i' % (case_color[i])
                if S_4_m[i] != 2. * h: print 'error on S_1_m'
                if S_1_m[i] != h: print 'error on S_2_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 13:
                print 'case_c=%i' % (case_color[i])
                if S_1_m[i] != h: print 'error on S_1_m=%f and h=%f' % (S_1_m[i], h)
                if S_4_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_3_m[i] != S_1_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i + 1])), 20)
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 14:
                print 'case_c=%i' % (case_color[i])
                if S_4_m[i] != h: print 'error on S_1_m'
                if S_3_m[i] != S_1_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i + 1])), 20)
                if S_1_m[i] != S_3_m[i - 1]: print 'error on S_2_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i - 1])), 20)

            if case_color[i] == 15:
                print 'case_c=%i' % (case_color[i])
                if S_4_m[i] != h: print 'error on S_1_m'
                if S_3_m[i] != S_1_m[i + 1]: print 'error on S_3_m: %.20f' % round((abs(S_3_m[i] - S_1_m[i + 1])), 20)
                if S_1_m[i] != S_3_m[i - 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i - 1])), 20)
            if case_color[i] == 16:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m=%f and h=%f' % (S_3_m[i], h)
                if S_4_m[i] != h: print 'error on S_4_m=%f and h=%f' % (S_4_m[i], h)
                if S_1_m[i] != S_3_m[i - 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i - 1])), 20)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
            if case_color[i] == 17:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_4_m[i] != 2. * h: print 'error on S_4_m=%.20f and h=%.20f' % (S_4_m[i], h)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 18:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_4_m[i] != h: print 'error on S_4_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 19:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != 2. * h: print 'error on S_3_m'
                if S_4_m[i] != h: print 'error on S_4_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 20:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m=%f and h=%f' % (S_3_m[i], h)
                if S_4_m[i] != h: print 'error on S_4_m=%f and h=%f' % (S_4_m[i], h)
                if S_2_m[i] != S_4_m[i + 1]: print 'error on S_2_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i + 1])), 20)
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 21:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_2_m[i] != S_4_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i + 1])), 20)
                if S_4_m[i] != S_2_m[i - 1]: print 'error on S_2_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i - 1])), 20)

            if case_color[i] == 22:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_2_m[i] != S_4_m[i + 1]: print 'error on S_4_m: %.20f' % round((abs(S_2_m[i] - S_4_m[i + 1])), 20)
                if S_4_m[i] != S_2_m[i - 1]: print 'error on S_2_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i - 1])), 20)
            if case_color[i] == 23:
                print 'case_c=%i' % (case_color[i])
                if S_2_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_3_m[i] != h: print 'error on S_3_m=%f and h=%f' % (S_3_m[i], h)
                if S_4_m[i] != S_2_m[i - 1]: print 'error on S_4_m: %.20f' % round((abs(S_4_m[i] - S_2_m[i - 1])), 20)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
            if case_color[i] == 24:
                print 'case_c=%i' % (case_color[i])
                if S_2_m[i] != h: print 'error on S_2_m'
                if S_3_m[i] != 2. * h: print 'error on S_3_m=%.20f and h=%.20f' % (S_3_m[i], h)
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 25:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_2_m[i] != h: print 'error on S_2_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 26:
                print 'case_c=%i' % (case_color[i])
                if S_2_m[i] != 2. * h: print 'error on S_2_m'
                if S_3_m[i] != h: print 'error on S_3_m'
                if S_6_m[i] != S_5_m[i + 1]: print 'error on S_6_m'
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if case_color[i] == 27:
                print 'case_c=%i' % (case_color[i])
                if S_3_m[i] != h: print 'error on S_3_m=%f and h=%f' % (S_3_m[i], h)
                if S_2_m[i] != h: print 'error on S_2_m=%f and h=%f' % (S_2_m[i], h)
                if S_1_m[i] != S_3_m[i + 1]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[i + 1])), 20)
                if S_5_m[i] != S_6_m[i - 1]: print 'error on S_5_m'
            if i == len(case_color) - 1:
                j = 0
            else:
                j = i + 1
            if case_color[i] == 28:
                print 'case_c=%i' % (case_color[i])
                if S_2_m[i] != h: print 'error on S_2_m'
                if S_1_m[i] != S_3_m[j]: print 'error on S_1_m: %.20f' % round((abs(S_1_m[i] - S_3_m[j])), 20)
                if S_3_m[i] != S_1_m[i - 1]: print 'error on S_3_m: %.20f' % round((abs(S_3_m[i] - S_1_m[i - 1])), 20)

        for i in range(len(S_5_m)):
            if S_6_m[i] != 0:
                if S_6_m[i] != S_5_m[i + 1]:
                    print 'error on S_6_m'
                    print 'case_c=%i' % (case_color[i])
                    print S_6_m[i], S_5_m[i + 1]
                else:
                    print 'S_6_m matched with S_5_m'
            if S_5_m[i] != 0:
                if S_5_m[i] != S_6_m[i - 1]:
                    print 'error on S_5_m'
                    print 'case_c=%i' % (case_color[i])
                    print S_6_m[i], S_5_m[i - 1]
                else:
                    print 'S_5_m matched with S_6_m'

    return h, radius, \
           x_inter, y_inter, \
           x_halcell, y_halcell, \
           x_cutcell, y_cutcell, \
           x_cc_m, y_cc_m, \
           S_1_m, S_2_m, S_3_m, S_4_m, S_5_m, S_6_m, \
           S_b_m, V_b_m, \
           case_color, \
           x_cc_tri, y_cc_tri, \
           S_1_tri, S_2_tri, S_3_tri, S_4_tri


if __name__ == '__main__':
    h_step = 0.0550
    # [   0.0213, 0.0251, 0.0274, 0.0355, 0.0405, 0.0430,
    #     0.0475, 0.0495, 0.0550, 0.0580, 0.0670, 0.0695]

    a = .3010 / 1.
    print a

    h, radius, \
    x_inter, y_inter, \
    x_halcell, y_halcell, \
    x_cutcell, y_cutcell, \
    x_cc_m, y_cc_m, \
    S_1_m, S_2_m, S_3_m, S_4_m, S_5_m, S_6_m, \
    S_b_m, V_b_m, \
    case_color, \
    x_cc_tri, y_cc_tri, \
    S_1_tri, S_2_tri, S_3_tri, S_4_tri = merging_cutcells(a, h_step, 301 * h_step, 301 * h_step)
