# -*- coding: utf-8 -*-
##
# \file     plot_boundary_cells.py
# \title    Shows and saves a map of the boundary cells (edges or corners).
# \author   Pierre Chobeau
# \version  0.1
# \date     2017, 28 Sep.
##
import numpy as np
from matplotlib import pyplot as plt

def plot_cells_circle(dx, bc_ed, bc_co, bc_in, ray, Nx, Ny):
    """
      Plot the boundary cells of the circular obstacle to check the corners
            and edges sent in the updates. Required to check any
            orverlaping or misplaced boundary cells during the calculation of the updates.

    :param  dx      spatial step after discretization, float (m).
    :param  bc_ed   array of edge cells that have 1 branch connected to the inner cell, 2D array of integers.
    :param  bc_co   array of corner cells that have 2 branches connected to the inner cells, 2D array of integers.
    :param  bc_in   array of inner cells that are zeroed out, 2D array of integers.
    :param  ray     radius of the circular obstacle, float (m).
    :param  Nx      discrete length of the domain (number of node) following the x-direction.
    :param  Ny      discrete length of the domain (number of node) following the y-direction.

    :return     Depiction of the boundary cells that map the edges and corners, and the inner cells.
    """
    fig = plt.figure('Boundary cells maps')
    from matplotlib.colors import from_levels_and_colors
    ax = fig.add_subplot(111)
    cmap, norm = from_levels_and_colors([0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5],
                                        [   'g' , 'y' , 'r' , 'b' , 'm' , 'c' ,'0.60','r' ])
    x_st_idx = int(round(Nx / 2. - (np.around(ray/dx) + 1)))
    x_en_idx = int(round(Nx / 2. + np.around(ray/dx) + 3))
    y_st_idx = int(round(Ny / 2. - (np.around(ray/dx) + 1)))
    y_en_idx = int(round(Ny / 2. + np.around(ray/dx) + 3))
    x = np.arange(x_st_idx,x_en_idx)
    y = np.arange(y_st_idx,y_en_idx)
    X, Y = np.meshgrid(x, y)
    plt.pcolormesh(X, Y, np.transpose(bc_co[x_st_idx:x_en_idx, y_st_idx:y_en_idx]), cmap=cmap, norm=norm)
    plt.pcolormesh(X, Y, np.transpose(bc_ed[x_st_idx:x_en_idx, y_st_idx:y_en_idx]), cmap=cmap, norm=norm)
    plt.pcolormesh(X, Y, np.transpose(bc_in[x_st_idx:x_en_idx, y_st_idx:y_en_idx]), cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(Nx / 2 - (np.around(ray/dx) + 1),Nx / 2 + np.around(ray/dx) + 3, 5))
    ax.set_yticks(np.arange(Ny / 2 - (np.around(ray/dx) + 1),Ny / 2 + np.around(ray/dx) + 3, 5))

    ax.grid(which='major', alpha=1.0)
    plt.grid()
    plt.xlim(Nx/2-(np.around(ray/dx)+1),Nx/2+np.around(ray/dx)+4)
    plt.ylim(Ny/2-(np.around(ray/dx)+1),Ny/2+np.around(ray/dx)+4)
    plt.title('bc_co ; bc_ed ; bc_in')
    plt.legend(('Corners', 'Edges ', 'Inside'), loc=1, fontsize=12)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def plot_cells_slope(bc_ed, bc_co, bc_in, slope_start, Nx, Ny):
    """
      Plot the boundary cells of the upward sloping part to check the corners
            and edges sent in the updates. Required to check any
            orverlaping or misplaced boundary cells during the calculation of the updates.

    :param  bc_ed   array of edge cells that have 1 branch connected to the inner cell, 2D array of integers.
    :param  bc_co   array of corner cells that have 2 branches connected to the inner cells, 2D array of integers.
    :param  bc_in   array of inner cells that are zeroed out, 2D array of integers.
    :param  slope_start     grid index at which the slope starts, integer.
    :param  Nx      discrete length of the domain (number of node) following the x-direction.
    :param  Ny      discrete length of the domain (number of node) following the y-direction.
    :return     Two plots of the boundary cells that map the edges and corners, and the inner cells.
    """
    fig = plt.figure('Boundary cells maps')
    from matplotlib.colors import from_levels_and_colors
    ax = fig.add_subplot(111)
    cmap, norm = from_levels_and_colors([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5],
                                        [  'g', 'y', 'r', 'b', 'm', 'c', '0.60', 'r'])
    x = np.arange(slope_start - 5, Nx)
    y = np.arange(0, Ny)
    X, Y = np.meshgrid(x, y)
    plt.pcolormesh(X, Y, np.transpose(bc_co[slope_start - 5: Nx, 0: Ny]), cmap=cmap, norm=norm)
    plt.pcolormesh(X, Y, np.transpose(bc_ed[slope_start - 5: Nx, 0: Ny]), cmap=cmap, norm=norm)
    plt.pcolormesh(X, Y, np.transpose(bc_in[slope_start - 5: Nx, 0: Ny]), cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(slope_start - 5, Nx, 10))
    ax.set_yticks(np.arange(0, Ny, 10))
    ax.grid(which='major', alpha=10.0)
    plt.grid()
    plt.xlim(slope_start - 2, Nx)
    plt.ylim(0, Ny)
    plt.title('bc_co ; bc_ed ; bc_in')
    plt.legend(('Corners', 'Edges', 'Inside'), loc=1, fontsize=12)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()