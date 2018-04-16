# -*- coding: utf-8 -*-
##
# \file     D_circ_final.py
# \title    Sort all the parameters out that are used in the updates
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 07 Sep.
##

import numpy as np
import C_circ_merge_cutcells as C_5

def D_circ_surf(h, radius, Lx, Ly):
    """
          Step 4/4, calculation of the cut-cells parameters - after merging - for the FVTD update.

    :param  radius  radius of the circle, float (m).
    :param  h       spatial step set, float (m).
    :param  Lx      length of the domain following the x axis, float (m).
    :param  Ly      length of the domain following the y axis, float (m).

    :param      x_halcell   x coordinates of the halved cells, list of floats (m).
    :param      x_halcell   y coordinates of the halved cells, list of floats (m).
    :param      x_cutcell   x coordinates of the cut cells, list of floats (m).
    :param      y_cutcell   y coordinates of the cut cells, list of floats (m).
    :param      x_cc   x coordinates of the new center of mass for the cut cells, list of floats (m).
    :param      y_cc   y coordinates of the new center of mass for the cut cells, list of floats (m).

    :return     the new parameters (volumes, surfaces, corrdinates) of the cut cells after merging.
    """
    dl = h
    Nx = np.int(np.round(Lx / dl));
    Ny = np.int(np.round(Ly / dl));
    ax = (np.floor(Nx / 2) + h) * dl;
    by = (np.floor(Ny / 2)) * dl;

    h, radius, \
    x_inter, y_inter, \
    x_halcell, y_halcell, \
    x_cutcell, y_cutcell, \
    x_cc, y_cc, \
    S_1, S_2, S_3, S_4, S_5, S_6, \
    S, V, \
    case_color, \
    x_cc_tri, y_cc_tri, \
    S_1_tri, S_2_tri, S_3_tri, S_4_tri = C_5.merging_cutcells(radius, h, ax, by)

    # ==============================================================================
    #   Nb == all solid cells for the   Staircases (sc),
    #                                   Halved-Corners (hc),
    #                                   Cut-Cells (cc).
    # ==============================================================================

    Nb_x = [];
    Nb_y = [];  # INNER-CELLS (p=0)
    #   DETECTION OF THE INNER-CELLS (SOLID)
    for i in range(int(ax / h - radius / h), int(ax / h + radius / h) + 1):
        for j in range(int(by / h - radius / h), int(by / h + radius / h) + 1):
            if (abs(i * h - ax) ** 2 + abs(j * h - by) ** 2) < (radius + 1 * h / 2) ** 2:
                Nb_x.append(i * h)
                Nb_y.append(j * h)

    Nb_x_idx = [int(round(x / h)) for x in Nb_x]
    Nb_y_idx = [int(round(y / h)) for y in Nb_y]

    #   CUT-CELLS COORDINATES BEFORE MERGING
    x_cutcell_idx = [int(round(x / h, 0)) for x in x_cutcell]
    y_cutcell_idx = [int(round(y / h, 0)) for y in y_cutcell]

    #   CUT-CELLS COORDINATES AFTER MERGING
    x_cc_idx = [int(round(x / h, 0)) for x in x_cc]
    y_cc_idx = [int(round(y / h, 0)) for y in y_cc]

    #   TRIANGULAR CELLS COORDINATES AFTER MERGING
    x_tri_idx = [int(round(x / h, 0)) for x in x_cc_tri]
    y_tri_idx = [int(round(y / h, 0)) for y in y_cc_tri]

    Nb_x_final_cc = []  # INNER-CELLS (p=0): for CUT-CELLS
    Nb_y_final_cc = []
    Nb_x_final_cc_idx = []
    Nb_y_final_cc_idx = []
    #   ELIMINATION OF THE OVERLAPPING INNER-CELLS (SOLID) & CUT-CELLS (USING indexes:idx)
    for i in range(len(Nb_x)):
        count = 0
        for j in range(len(x_cc_idx)):
            if Nb_x_idx[i] == x_cc_idx[j] and Nb_y_idx[i] == y_cc_idx[j]:
                count += 1
        if count == 0:
            Nb_x_final_cc.append(Nb_x[i])
            Nb_y_final_cc.append(Nb_y[i])
            Nb_x_final_cc_idx.append(Nb_x_idx[i])
            Nb_y_final_cc_idx.append(Nb_y_idx[i])

    x_halcell_idx = [int(round(x / h, 0)) for x in x_halcell]
    y_halcell_idx = [int(round(y / h, 0)) for y in y_halcell]
    Nb_x_final_hc = []  # INNER-CELLS (p=0): for 1/2-CORNERS
    Nb_y_final_hc = []
    Nb_x_final_hc_idx = []
    Nb_y_final_hc_idx = []
    #   ELIMINATION OF THE OVERLAPPING INNER-CELLS (SOLID) & 1/2-CORNERS (USING indexes:idx)
    for i in range(len(Nb_x)):
        count = 0
        for j in range(len(x_halcell)):
            if Nb_x_idx[i] == x_halcell_idx[j] and Nb_y_idx[i] == y_halcell_idx[j]:
                count += 1
        if count == 0:
            Nb_x_final_hc.append(Nb_x[i])
            Nb_y_final_hc.append(Nb_y[i])
            Nb_x_final_hc_idx.append(Nb_x_idx[i])
            Nb_y_final_hc_idx.append(Nb_y_idx[i])

    Nb_x_final_hc_idx = [int(round(x / h)) for x in Nb_x_final_hc]
    Nb_y_final_hc_idx = [int(round(y / h)) for y in Nb_y_final_hc]
    x_hal_corners = [];
    x_hal_corners_idx = []
    y_hal_corners = [];
    y_hal_corners_idx = []
    x_hal_edges = [];
    x_hal_edges_idx = []
    y_hal_edges = [];
    y_hal_edges_idx = []
    x_corners_idx_1st = []
    y_corners_idx_1st = []
    x_corners_idx_2nd = []
    y_corners_idx_2nd = []
    x_corners_idx_3rd = []
    y_corners_idx_3rd = []
    x_corners_idx_4th = []
    y_corners_idx_4th = []
    x_edges_idx_1st = []
    y_edges_idx_1st = []
    x_edges_idx_2nd = []
    y_edges_idx_2nd = []
    x_edges_idx_3rd = []
    y_edges_idx_3rd = []
    x_edges_idx_4th = []
    y_edges_idx_4th = []
    #   DETECTION OF THE CORNERS for the 1/2-corners (green stars)
    for i in range(len(x_halcell_idx)):
        count = 0
        for j in range(len(Nb_x_final_hc_idx)):
            if (x_halcell_idx[i] + 1) == Nb_x_final_hc_idx[j] and y_halcell_idx[i] == Nb_y_final_hc_idx[j]:
                count += 1
            if (x_halcell_idx[i] - 1) == Nb_x_final_hc_idx[j] and y_halcell_idx[i] == Nb_y_final_hc_idx[j]:
                count += 1
            if x_halcell_idx[i] == Nb_x_final_hc_idx[j] and (y_halcell_idx[i] + 1) == Nb_y_final_hc_idx[j]:
                count += 1
            if x_halcell_idx[i] == Nb_x_final_hc_idx[j] and (y_halcell_idx[i] - 1) == Nb_y_final_hc_idx[j]:
                count += 1

        if count == 2:
            x_hal_corners.append(x_halcell[i])
            y_hal_corners.append(y_halcell[i])
            x_hal_corners_idx.append(x_halcell_idx[i])
            y_hal_corners_idx.append(y_halcell_idx[i])
        else:
            x_hal_edges.append(x_halcell[i])
            y_hal_edges.append(y_halcell[i])
            x_hal_edges_idx.append(x_halcell_idx[i])
            y_hal_edges_idx.append(y_halcell_idx[i])

    x_S34_idx = [];
    y_S34_idx = [];
    x_S23_idx = [];
    y_S23_idx = []
    x_S12_idx = [];
    y_S12_idx = [];
    x_S41_idx = [];
    y_S41_idx = []
    x_S1_idx = [];
    y_S1_idx = [];
    x_S2_idx = [];
    y_S2_idx = []
    x_S3_idx = [];
    y_S3_idx = [];
    x_S4_idx = [];
    y_S4_idx = []
    # ==============================================================================
    #   Gives the coordinates of the cells adjacent to the boundary cells
    #   The sorting distinguisches between:
    #                - the 1-cell connections (x_S_3, x_S_4 ...) and
    #                - the 2-cell connections (x_S_34, x_S_23 ... ).
    # ==============================================================================
    x_c_1_4 = [x_hal_corners_idx[i] for i in range(len(x_hal_corners_idx) / 4)]
    y_c_1_4 = [y_hal_corners_idx[i] for i in range(len(y_hal_corners_idx) / 4)]
    x_c_2_4 = [x_hal_corners_idx[i] for i in range(len(x_hal_corners_idx) / 4, len(x_hal_corners_idx) / 2)]
    y_c_2_4 = [y_hal_corners_idx[i] for i in range(len(y_hal_corners_idx) / 4, len(x_hal_corners_idx) / 2)]
    x_c_3_4 = [x_hal_corners_idx[i] for i in range(len(x_hal_corners_idx) / 2, 3 * len(x_hal_corners_idx) / 4)]
    y_c_3_4 = [y_hal_corners_idx[i] for i in range(len(y_hal_corners_idx) / 2, 3 * len(x_hal_corners_idx) / 4)]
    x_c_4_4 = [x_hal_corners_idx[i] for i in range(3 * len(x_hal_corners_idx) / 4, len(x_hal_corners_idx))]
    y_c_4_4 = [y_hal_corners_idx[i] for i in range(3 * len(y_hal_corners_idx) / 4, len(x_hal_corners_idx))]

    i_1_4 = [];
    #    print x_c_1_4, y_c_1_4

    #   get the 1-cell connection coordinates and
    #   eliminate those from the list using i_1_4

    # SORT OUT THE EDGE CELLS ACCORDING TO THE QUADRANTS
    for idx in range(len(x_hal_edges_idx)):
        r_hypo = np.sqrt((x_hal_edges[idx]-ax)**2+(y_hal_edges[idx]-by)**2)
        if  (np.arcsin((x_hal_edges[idx] - ax)/r_hypo)) > np.pi/4. and\
            (np.arcsin((x_hal_edges[idx] - ax)/r_hypo)) < 3.*np.pi/4.:
            x_edges_idx_1st.append(x_hal_edges_idx[idx])
            y_edges_idx_1st.append(y_hal_edges_idx[idx])
        if  (np.arcsin((y_hal_edges[idx] - by)/r_hypo)) > np.pi/4. and\
            (np.arcsin((y_hal_edges[idx] - by)/r_hypo)) < 3.*np.pi/4.:
            x_edges_idx_2nd.append(x_hal_edges_idx[idx])
            y_edges_idx_2nd.append(y_hal_edges_idx[idx])
        if  (np.arcsin((ax - x_hal_edges[idx])/r_hypo)) > np.pi/4. and\
            (np.arcsin((ax - x_hal_edges[idx])/r_hypo)) < 3.*np.pi/4.:
            x_edges_idx_3rd.append(x_hal_edges_idx[idx])
            y_edges_idx_3rd.append(y_hal_edges_idx[idx])
        if  (np.arcsin((by - y_hal_edges[idx])/r_hypo)) > np.pi/4. and\
            (np.arcsin((by - y_hal_edges[idx])/r_hypo)) < 3.*np.pi/4.:
            x_edges_idx_4th.append(x_hal_edges_idx[idx])
            y_edges_idx_4th.append(y_hal_edges_idx[idx])

    # SORT OUT THE CORNER CELLS ACCORDING TO THE QUADRANTS
    if x_hal_corners:
        for idx in range(len(x_hal_corners_idx)):
            i = x_hal_corners_idx[idx]
            j = y_hal_corners_idx[idx]

            #           1st quadrant
            if ((y_hal_corners[idx] - by) > 0 and (x_hal_corners[idx] - ax) > 0):
                for jdx in range(len(x_hal_edges)):
                    x_corners_idx_1st.append(x_hal_corners_idx[idx])
                    y_corners_idx_1st.append(y_hal_corners_idx[idx])
                    if (
                                ((x_hal_corners_idx[idx] - 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] + 1) == y_hal_edges_idx[jdx])):
                        x_S4_idx.append(i)
                        y_S4_idx.append(j + 1)
                        i_1_4.append(idx)

                    elif (
                                ((x_hal_corners_idx[idx] + 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] - 1) == y_hal_edges_idx[jdx])):
                        x_S1_idx.append(i + 1)
                        y_S1_idx.append(j)
                        i_1_4.append(idx)

                        #           2nd quadrant
            if ((y_hal_corners[idx] - by) > 0 and (x_hal_corners[idx] - ax) < 0):
                for jdx in range(len(x_hal_edges)):
                    x_corners_idx_2nd.append(x_hal_corners_idx[idx])
                    y_corners_idx_2nd.append(y_hal_corners_idx[idx])
                    if (
                                ((x_hal_corners_idx[idx] + 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] + 1) == y_hal_edges_idx[jdx])):
                        x_S4_idx.append(i)
                        y_S4_idx.append(j + 1)
                    # i_1_4.append(i)

                    elif (
                                ((x_hal_corners_idx[idx] - 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] - 1) == y_hal_edges_idx[jdx])):
                        x_S3_idx.append(i - 1)
                        y_S3_idx.append(j)
                        #                            i_1_4.append(i)

                        #           3rd quadrant
            if ((y_hal_corners[idx] - by) < 0 and (x_hal_corners[idx] - ax) < 0):
                for jdx in range(len(x_hal_edges)):
                    x_corners_idx_3rd.append(x_hal_corners_idx[idx])
                    y_corners_idx_3rd.append(y_hal_corners_idx[idx])
                    if (
                                ((x_hal_corners_idx[idx] + 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] - 1) == y_hal_edges_idx[jdx])):
                        x_S2_idx.append(i)
                        y_S2_idx.append(j - 1)
                    # i_1_4.append(i)

                    elif (
                                ((x_hal_corners_idx[idx] - 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] + 1) == y_hal_edges_idx[jdx])):
                        x_S3_idx.append(i - 1)
                        y_S3_idx.append(j)
                        #                            i_1_4.append(i)

                        #           4th quadrant
            if ((y_hal_corners[idx] - by) < 0 and (x_hal_corners[idx] - ax) > 0):
                for jdx in range(len(x_hal_edges)):
                    x_corners_idx_4th.append(x_hal_corners_idx[idx])
                    y_corners_idx_4th.append(y_hal_corners_idx[idx])
                    if (
                                ((x_hal_corners_idx[idx] - 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] - 1) == y_hal_edges_idx[jdx])):
                        x_S2_idx.append(i)
                        y_S2_idx.append(j - 1)
                    # i_1_4.append(i)

                    elif (
                                ((x_hal_corners_idx[idx] + 1) == x_hal_edges_idx[jdx]) and
                                ((y_hal_corners_idx[idx] + 1) == y_hal_edges_idx[jdx])):
                        x_S1_idx.append(i + 1)
                        y_S1_idx.append(j)

                        #   get the 2-cell connections:
                        #               1/. for single 2-cell ceonnection cases (first if ...),
                        #       or      2/. for multiple 2-cell connection cases (after the else: ...).
                        #        print len(x_c_1_4)
        if len(x_c_1_4) > 1:
            if len(x_c_1_4) == 2:
                x_c_1_4 = x_c_1_4[0]
                y_c_1_4 = y_c_1_4[0]
                x_S41_idx.append(x_c_1_4)
                y_S41_idx.append(y_c_1_4 + 1)

                x_c_2_4 = x_c_2_4[0]
                y_c_2_4 = y_c_2_4[0]
                x_S34_idx.append(x_c_2_4 - 1)
                y_S34_idx.append(y_c_2_4)

                x_c_3_4 = x_c_3_4[0]
                y_c_3_4 = y_c_3_4[0]
                x_S23_idx.append(x_c_3_4)
                y_S23_idx.append(y_c_3_4 - 1)

                x_c_4_4 = x_c_4_4[0]
                y_c_4_4 = y_c_4_4[0]
                x_S12_idx.append(x_c_4_4 + 1)
                y_S12_idx.append(y_c_4_4)
            else:
                x_c_1_4[:] = (x_c_1_4[i] for i in range(len(x_c_1_4)) if not i in i_1_4)
                y_c_1_4[:] = (y_c_1_4[i] for i in range(len(y_c_1_4)) if not i in i_1_4)
                x_c_2_4[:] = (x_c_2_4[i] for i in range(len(x_c_2_4)) if not i in i_1_4)
                y_c_2_4[:] = (y_c_2_4[i] for i in range(len(y_c_2_4)) if not i in i_1_4)
                x_c_3_4[:] = (x_c_3_4[i] for i in range(len(x_c_3_4)) if not i in i_1_4)
                y_c_3_4[:] = (y_c_3_4[i] for i in range(len(y_c_3_4)) if not i in i_1_4)
                x_c_4_4[:] = (x_c_4_4[i] for i in range(len(x_c_4_4)) if not i in i_1_4)
                y_c_4_4[:] = (y_c_4_4[i] for i in range(len(y_c_4_4)) if not i in i_1_4)

                for i in range(0, len(x_c_1_4)):
                    x_S41_idx.append(x_c_1_4[i] + 1)
                    y_S41_idx.append(y_c_1_4[i])
                    x_S41_idx.append(x_c_1_4[i])
                    y_S41_idx.append(y_c_1_4[i] + 1)

                    x_S34_idx.append(x_c_2_4[i])
                    y_S34_idx.append(y_c_2_4[i] + 1)
                    x_S34_idx.append(x_c_2_4[i] - 1)
                    y_S34_idx.append(y_c_2_4[i])

                    x_S23_idx.append(x_c_3_4[i] - 1)
                    y_S23_idx.append(y_c_3_4[i])
                    x_S23_idx.append(x_c_3_4[i])
                    y_S23_idx.append(y_c_3_4[i] - 1)

                    x_S12_idx.append(x_c_4_4[i])
                    y_S12_idx.append(y_c_4_4[i] - 1)
                    x_S12_idx.append(x_c_4_4[i] + 1)
                    y_S12_idx.append(y_c_4_4[i])

    # ==============================================================================
    #   Set the Boundary-Cell surfaces for the 1/2-corners
    # ==============================================================================
    S1_hc = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    S2_hc = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    S3_hc = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    S4_hc = np.zeros((Nx + 1, Ny + 1), dtype=np.float64)
    if x_hal_corners:
        for idx in range(len(x_hal_corners_idx)):
            i = x_hal_corners_idx[idx]
            j = y_hal_corners_idx[idx]
            #           1st quadrant
            if ((y_hal_corners[idx] - by) > 0 and (x_hal_corners[idx] - ax) > 0):
                S2_hc[i, j] = h
                S3_hc[i, j] = h
            # 2nd quadrant
            if ((y_hal_corners[idx] - by) > 0 and (x_hal_corners[idx] - ax) < 0):
                S1_hc[i, j] = h
                S2_hc[i, j] = h
            # 3rd quadrant
            if ((y_hal_corners[idx] - by) < 0 and (x_hal_corners[idx] - ax) < 0):
                S1_hc[i, j] = h
                S4_hc[i, j] = h
            # 4th quadrant
            if ((y_hal_corners[idx] - by) < 0 and (x_hal_corners[idx] - ax) > 0):
                S3_hc[i, j] = h
                S4_hc[i, j] = h

    # ==============================================================================
    #   CALCULATION OF THE CIRCUMFERENCES & INNER-SURFACES
    # ==============================================================================
    Surf_in_sc = len(Nb_x_final_hc) * h ** 2
    Surf_in_hc = len(Nb_x_final_hc) * h ** 2 + len(x_hal_corners) * (h ** 2) / 2.
    Surf_in_cc = len(Nb_x_final_cc) * h ** 2 + len(x_cutcell) * h ** 2 - sum(V)
    Surf_in_an = np.pi * radius ** 2
    Surf_in_sc_rate = Surf_in_an / Surf_in_sc
    Surf_in_hc_rate = Surf_in_an / Surf_in_hc
    Surf_in_cc_rate = Surf_in_an / Surf_in_cc

    Circ_in_sc = len(x_hal_edges) * h + len(x_hal_corners) * 2 * h
    Circ_in_hc = len(x_hal_edges) * h + len(x_hal_corners) * h * np.sqrt(2)
    Circ_in_cc = sum(S)
    Circ_in_an = 2 * np.pi * radius
    Circ_in_sc_rate = Circ_in_an / Circ_in_sc
    Circ_in_hc_rate = Circ_in_an / Circ_in_hc
    Circ_in_cc_rate = Circ_in_an / Circ_in_cc

    phi = np.linspace(0, 2 * np.pi, 10000)
    x_circ = np.zeros((len(phi)), dtype=np.float64)
    x_circ = radius * np.cos(phi) + ax
    y_circ = np.zeros((len(phi)), dtype=np.float64)
    y_circ = radius * np.sin(phi) + by
    #
    import matplotlib.pyplot as plt
    plot_1 = False
    if plot_1:
        ##==============================================================================
        ##   CHECK THE EDGES FOR THE TLM 1ST - 4TH...
        ##==============================================================================
        fig = plt.figure('TLM edges and corners')
        ax1 = fig.add_subplot(111)
        ax1.set_aspect('equal')

        plt.plot(ax, by, 'k+', markersize=10, markeredgewidth=4, markeredgecolor='k', markerfacecolor='None')

        plt.plot(x_circ, y_circ, 'r-')

        plt.plot([h for i in range(Ny)], [i * h for i in range(Ny)], 'r*', markersize=10,
                 markeredgewidth=4,
                 markeredgecolor='r', markerfacecolor='None')

        plt.plot([i*h for i in x_edges_idx_1st], [i*h for i in y_edges_idx_1st], 'cs', markersize=10, markeredgewidth=4,
                 markeredgecolor='c', markerfacecolor='None')

        plt.plot([i*h for i in x_edges_idx_2nd], [i*h for i in y_edges_idx_2nd], 'ro', markersize=8, markeredgewidth=2,
                 markeredgecolor='r', markerfacecolor='None')

        plt.plot([i*h for i in x_edges_idx_3rd], [i*h for i in y_edges_idx_3rd], 'g*', markersize=20, markeredgewidth=2,
                 markeredgecolor='g', markerfacecolor='None')

        plt.plot([i*h for i in x_edges_idx_4th], [i*h for i in y_edges_idx_4th], 'yd', markersize=20, markeredgewidth=2,
                 markeredgecolor='y', markerfacecolor='None')

        plt.plot([i*h for i in x_corners_idx_1st], [i*h for i in y_corners_idx_1st], 'cs', markersize=2, markeredgewidth=2,
                 markeredgecolor='c', markerfacecolor='None')

        plt.plot([i*h for i in x_corners_idx_2nd], [i*h for i in y_corners_idx_2nd], 'ro', markersize=2, markeredgewidth=2,
                 markeredgecolor='r', markerfacecolor='None')

        plt.plot([i*h for i in x_corners_idx_3rd], [i*h for i in y_corners_idx_3rd], 'g*', markersize=2, markeredgewidth=2,
                 markeredgecolor='g', markerfacecolor='None')

        plt.plot([i*h for i in x_corners_idx_4th], [i*h for i in y_corners_idx_4th], 'yd', markersize=2, markeredgewidth=2,
                 markeredgecolor='y', markerfacecolor='None')

        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=3.)
        # plt.xlim(ax - 1.2 * radius, ax + 1.4 * radius)
        # plt.ylim(by - 1.2 * radius, by + 1.4 * radius)

    ##==============================================================================
    ##   Figure for the Cut-Cell (cc) approach
    ##==============================================================================
        fig = plt.figure('GRID vs. CIRCLE (3): CUT-CELL')
        ax1 = fig.add_subplot(111)
        ax1.set_aspect('equal')

        plt.plot(ax, by, 'k+', markersize=10, markeredgewidth=4, markeredgecolor='k', markerfacecolor='None')

        plt.plot(x_circ, y_circ, 'r-')

        plt.plot(Nb_x_final_cc, Nb_y_final_cc, 'cs', markersize=10, markeredgewidth=4,
                 markeredgecolor='c', markerfacecolor='None')

        plt.plot(x_cutcell, y_cutcell, 'ro', markersize=8, markeredgewidth=2,
                 markeredgecolor='r', markerfacecolor='None')

        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=3.)
        plt.xlim(ax - 1.2 * radius, ax + 1.4 * radius)
        plt.ylim(by - 1.2 * radius, by + 1.4 * radius)

        plt.show()

    # ==============================================================================
    #   Figure for the Halved-Cell (hc) approach
    # ==============================================================================
    plot_2 = False
    if plot_2:
        fig = plt.figure('GRID vs. CIRCLE (4): 1/2-CELL, h=%f' % (h))
        ax1 = fig.add_subplot(111)
        ax1.set_aspect('equal')

        plt.plot(ax, by, 'k+', markersize=10, markeredgewidth=4, markeredgecolor='k', markerfacecolor='None')

        plt.plot(x_circ, y_circ, 'r-')

        plt.plot(Nb_x_final_hc, Nb_y_final_hc, 'cs', markersize=10, markeredgewidth=4,
                 markeredgecolor='c', markerfacecolor='None')

        plt.plot(x_halcell, y_halcell, 'ro', markersize=8, markeredgewidth=2,
                 markeredgecolor='r', markerfacecolor='None')

        plt.plot(x_hal_corners, y_hal_corners, 'g*', markersize=20, markeredgewidth=2,
                 markeredgecolor='g', markerfacecolor='None')

        plt.plot(x_hal_edges, y_hal_edges, 'yd', markersize=20, markeredgewidth=2,
                 markeredgecolor='y', markerfacecolor='None')

        x_S4_idx_h = [i * h for i in x_S4_idx]
        y_S4_idx_h = [i * h for i in y_S4_idx]
        plt.plot(x_S4_idx_h, y_S4_idx_h, 'g+', markersize=10, markeredgewidth=3,
                 markeredgecolor='g', markerfacecolor='None')
        x_S3_idx_h = [i * h for i in x_S3_idx]
        y_S3_idx_h = [i * h for i in y_S3_idx]
        plt.plot(x_S3_idx_h, y_S3_idx_h, 'bx', markersize=10, markeredgewidth=3,
                 markeredgecolor='b', markerfacecolor='None')
        x_S2_idx_h = [i * h for i in x_S2_idx]
        y_S2_idx_h = [i * h for i in y_S2_idx]
        plt.plot(x_S2_idx_h, y_S2_idx_h, 'm+', markersize=10, markeredgewidth=3,
                 markeredgecolor='m', markerfacecolor='None')
        x_S1_idx_h = [i * h for i in x_S1_idx]
        y_S1_idx_h = [i * h for i in y_S1_idx]
        plt.plot(x_S1_idx_h, y_S1_idx_h, 'rx', markersize=10, markeredgewidth=3,
                 markeredgecolor='r', markerfacecolor='None')

        x_S12_idx_h = [i * h for i in x_S12_idx]
        y_S12_idx_h = [i * h for i in y_S12_idx]
        plt.plot(x_S12_idx_h, y_S12_idx_h, 'r*', markersize=10, markeredgewidth=2,
                 markeredgecolor='r', markerfacecolor='None')
        x_S23_idx_h = [i * h for i in x_S23_idx]
        y_S23_idx_h = [i * h for i in y_S23_idx]
        plt.plot(x_S23_idx_h, y_S23_idx_h, 'k*', markersize=10, markeredgewidth=2,
                 markeredgecolor='k', markerfacecolor='None')
        x_S34_idx_h = [i * h for i in x_S34_idx]
        y_S34_idx_h = [i * h for i in y_S34_idx]
        plt.plot(x_S34_idx_h, y_S34_idx_h, 'b*', markersize=10, markeredgewidth=2,
                 markeredgecolor='b', markerfacecolor='None')
        x_S41_idx_h = [i * h for i in x_S41_idx]
        y_S41_idx_h = [i * h for i in y_S41_idx]
        plt.plot(x_S41_idx_h, y_S41_idx_h, 'g*', markersize=10, markeredgewidth=2,
                 markeredgecolor='g', markerfacecolor='None')
        plt.plot(x_inter, y_inter, 'ro--')

        ax1.set_xticks(np.arange(h / 2, 2 * ax + h + h / 2, h), minor=True)
        ax1.set_yticks(np.arange(h / 2, 2 * by + h + h / 2, h), minor=True)
        plt.grid(True, which='minor', color='k', linestyle='-', lw=3.)
        plt.xlim(ax - 1.2 * radius, ax + 1.4 * radius)
        plt.ylim(by - 1.2 * radius, by + 1.4 * radius)
        plt.show()
    # print 'Circle drawn'
    return  x_hal_corners_idx, y_hal_corners_idx,\
            x_hal_edges_idx, y_hal_edges_idx,\
            Nb_x_final_hc_idx, Nb_y_final_hc_idx,\
            x_corners_idx_1st, y_corners_idx_1st, x_corners_idx_2nd, y_corners_idx_2nd,\
            x_corners_idx_3rd, y_corners_idx_3rd, x_corners_idx_4th, y_corners_idx_4th, \
            x_edges_idx_1st, y_edges_idx_1st, x_edges_idx_2nd, y_edges_idx_2nd, \
            x_edges_idx_3rd, y_edges_idx_3rd, x_edges_idx_4th, y_edges_idx_4th


def plot_h_range(h, radius, Lx, Ly):
    Surf_SC = [];
    Surf_HC = [];
    Surf_CC = []
    Circ_SC = [];
    Circ_HC = [];
    Circ_CC = []

    for i in range(len(h)):
        Surf_in_sc_rate, Surf_in_hc_rate, Surf_in_cc_rate, \
        Circ_in_sc_rate, Circ_in_hc_rate, Circ_in_cc_rate, \
        x_halcell_idx, y_halcell_idx, \
        x_hal_corners_idx, y_hal_corners_idx, \
        x_hal_edges_idx, y_hal_edges_idx, \
        Nb_x_final_hc_idx, Nb_y_final_hc_idx, \
        x_cc_idx, y_cc_idx, \
        x_cutcell_idx, y_cutcell_idx, \
        Nb_x_final_cc_idx, Nb_y_final_cc_idx, \
        S_1, S_2, S_3, S_4, S_5, S_6, S, V, \
        case_color, \
        x_inter, y_inter, \
        x_S1_idx, y_S1_idx, x_S2_idx, y_S2_idx, \
        x_S3_idx, y_S3_idx, x_S4_idx, y_S4_idx, \
        x_S12_idx, y_S12_idx, x_S23_idx, y_S23_idx, \
        x_S34_idx, y_S34_idx, x_S41_idx, y_S41_idx, \
        S1_hc, S2_hc, S3_hc, S4_hc = D_circ_surf(h[i], radius, Lx, Ly)
        Surf_SC.append(Surf_in_sc_rate)
        Surf_HC.append(Surf_in_hc_rate)
        Surf_CC.append(Surf_in_cc_rate)
        Circ_SC.append(Circ_in_sc_rate)
        Circ_HC.append(Circ_in_hc_rate)
        Circ_CC.append(Circ_in_cc_rate)

    plot = False
    if plot:
        import matplotlib.pyplot as plt
        plt.figure('Relative Surface')
        plt.plot(h, Surf_SC, 'gs:', markersize=4, markeredgewidth=1, markeredgecolor='g', markerfacecolor='None')
        plt.plot(h, Surf_HC, 'r^--', markersize=4, markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
        plt.plot(h, Surf_CC, 'ko-', markersize=4, markeredgewidth=1, markeredgecolor='k', markerfacecolor='None')
        plt.grid()
        plt.legend((
            'Staircaises',
            '1/2-corners',
            'Cut-cells',
        ), fontsize=14, loc=2)

        plt.figure('Relative Circumference')
        plt.plot(h, Circ_SC, 'gs:', markersize=4, markeredgewidth=1, markeredgecolor='g', markerfacecolor='None')
        plt.plot(h, Circ_HC, 'r^--', markersize=4, markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
        plt.plot(h, Circ_CC, 'ko-', markersize=4, markeredgewidth=1, markeredgecolor='k', markerfacecolor='None')
        plt.grid()
        plt.legend((
            'Staircaises',
            '1/2-corners',
            'Cut-cells',
        ), fontsize=14, loc=4)
    diff_sc = [abs(x - 1) for x in Circ_SC]
    diff_hc = [abs(x - 1) for x in Circ_HC]
    diff_cc = [abs(x - 1) for x in Circ_CC]
    return diff_sc, diff_hc, diff_cc

if __name__ == '__main__':
    # diff_sc, diff_sc, diff_sc=\
    # plot_h_range([  0.0213, 0.0251, 0.0274, 0.0355, 0.0405, 0.0430,
    #                 0.0475, 0.0495, 0.0550, 0.0580, 0.0670, 0.0695
    #                ],.301,20.,20.)
    # h_range = [  0.0213, 0.0251, 0.0274, 0.0355, 0.0405, 0.0430,
    #                 0.0475, 0.0495, 0.0550, 0.0580, 0.0670, 0.0695 ]

    diff_sc, diff_sc, diff_sc = \
        plot_h_range([0.0405], .301, 20., 20.)
    h_range = [0.0405]

    for h in h_range:
        print h
        h_step = h
        Surf_in_sc_rate, Surf_in_hc_rate, Surf_in_cc_rate, \
        Circ_in_sc_rate, Circ_in_hc_rate, Circ_in_cc_rate, \
        x_halcell_idx, y_halcell_idx, \
        x_hal_corners_idx, y_hal_corners_idx, \
        x_hal_edges_idx, y_hal_edges_idx, \
        Nb_x_final_hc_idx, Nb_y_final_hc_idx, \
        x_cc_idx, y_cc_idx, \
        x_cutcell_idx, y_cutcell_idx, \
        Nb_x_final_cc_idx, Nb_y_final_cc_idx, \
        S_1, S_2, S_3, S_4, S_5, S_6, S, V, \
        case_color, \
        x_inter, y_inter, \
        x_S1_idx, y_S1_idx, x_S2_idx, y_S2_idx, \
        x_S3_idx, y_S3_idx, x_S4_idx, y_S4_idx, \
        x_S12_idx, y_S12_idx, x_S23_idx, y_S23_idx, \
        x_S34_idx, y_S34_idx, x_S41_idx, y_S41_idx, \
        S1_hc, S2_hc, S3_hc, S4_hc \
            = D_circ_surf(h_step, .301, 101 * h_step, 101 * h_step)