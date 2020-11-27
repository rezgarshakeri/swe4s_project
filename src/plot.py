"""This file includes function for plotting
the Temperature distribution in a 2d domain.
"""

import numpy as np
import matplotlib.pyplot as plt
import FE_subroutines as fe


def plot_temp(nelx, nely, T0_bottom, T0_left, K, F, grid):
    """This function plots the Temperature distribution
    in a 2d domain.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction.
    nely:       integer
                number of elements in the y direction.
    T0_bottom:  float
                prescribed temperature on the bottom edge
    T0_left:    float
                prescribed temperature on the left edge
    K:          float (2d array)
                assembled stiffness matrix
    F:          float (1d array)
                assembled forcing vector
    grid:       boolean
                true, if the user chooses to plot the grid

    Output:
    -------
    A contour plot for the temperature distriboution which
    looks like the following sketch for nelx = 2, nely = 2
    6---------7----------8
    |         |   (3)    |
    |   (2)   |      ----5
    |      ---4-----/    |
    3-----/   |   (1)    |
    |         |      ----2
    |   (0)   |     /
    |     ----1----/
    0----/
    """

    # Calls functions to return required output for the plots
    d0, ID, LM = fe.setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    d = fe.solvedr(nelx, nely, K, F, d0)
    nel, lpx, lpy, nnp, ndof, nen, neq = fe.setup(nelx, nely)
    x, y = fe.phys_coord(nelx, nely)
    IEN = fe.connectivity(nelx, nely)

    # Creates space for a figure to be drawn
    fig, ax = plt.subplots()

    # Arrange the temperature array in the physical ordering
    d1 = np.zeros(nnp)
    for i in range(nnp):
        d1[i] = d[ID[i][0]][0]

    # Creat and plot elements/grid
    if grid:
        xx = np.zeros(nen + 1)
        yy = np.zeros(nen + 1)
        dd = np.zeros(nen + 1)
        for i in range(nel):
            for j in range(nen):
                xx[j] = x[IEN[j, i]]
                yy[j] = y[IEN[j, i]]
                dd[j] = d1[IEN[j, i]]
            xx[nen] = x[IEN[0, i]]
            yy[nen] = y[IEN[0, i]]
            dd[nen] = d1[IEN[0, i]]
            plt.fill(xx, yy, edgecolor='black', fill=False)

    # Create and populate 2D Grid and the node's corresponding temperature
    X = np.zeros((lpx, lpy))
    Y = np.zeros((lpx, lpy))
    D = np.zeros((lpx, lpy))
    for i in range(lpx):
        for j in range(lpy):
            X[i][j] = x[j*lpx + i][0]
            Y[i][j] = y[j*lpx + i][0]
            D[i][j] = d1[j*lpx + i]

    # Turn off the right and top sides of the bounding box
    right_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side = ax.spines["right"]
    top_side.set_visible(False)

    # Plot
    plt.contourf(X, Y, D, cmap=plt.get_cmap('coolwarm'))
    plt.colorbar()
    plt.title('Temperature distribution')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    fig.tight_layout()
    plt.savefig('fig.png', bbox_inches='tight')

    return
