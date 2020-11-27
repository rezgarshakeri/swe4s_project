import matplotlib.pyplot as plt
import numpy as np

import FE_subroutines as fe 





def plot(nelx, nely, T0_bottom, T0_left, K, F):

    d0, ID, LM = fe.setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    d = fe.solvedr(nelx, nely, K, F, d0)
    nel, lpx, lpy, nnp, ndof, nen, neq = fe.setup(nelx, nely)
    x, y = fe.phys_coord(nelx, nely)
    IEN = fe.connectivity(nelx, nely)

    # creates space for a figure to be drawn 
    fig = plt.figure()
    
    d1 = np.zeros(nnp)
    for i in range(nnp):
        d1[i] = d[ID[i][0]][0]

    xx = np.zeros(nen + 1)
    yy = np.zeros(nen + 1)
    dd = np.zeros(nen + 1)
    for i in range(nel):
        for j in range(nen):
            xx[j] = x[ IEN[j, i] ]
            yy[j] = y[ IEN[j, i] ]
            dd[j] = d1[ IEN[j, i] ]
        xx[nen] = x[IEN[0, i]]
        yy[nen] = y[IEN[0, i]]
        dd[nen] = d1[IEN[0, i]]
        plt.fill(xx, yy, edgecolor='black', fill=False)

    X = np.zeros((lpx, lpy))
    Y = np.zeros((lpx, lpy))
    D = np.zeros((lpx, lpy))
    for i in range(lpx):
        for j in range(lpy):
            X[i][j] = x [ j*lpx + i ][0]
            Y[i][j] = y [ j*lpx + i ][0]
            D[i][j] = d1[ j*lpx + i ]

    plt.contourf(X, Y, D, cmap=plt.get_cmap('coolwarm'))
    plt.colorbar()

    plt.title('Temperature distribution')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    fig.tight_layout()
    plt.savefig('fig.png', bbox_inches='tight')


        
