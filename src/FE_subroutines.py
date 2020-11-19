""" This file includes the necessary modules for
implementing Finite Element Method.
"""

import numpy as np
import math
import sys


def setup(nelx, nely):
    # Total number of elements in the domain
    nel = nelx*nely

    # Number of nodes in the x direction
    lpx = nelx + 1

    # Number of nodes in the y direction
    lpy = nely + 1

    # Total number of nodes in the domain
    nnp = lpx*lpy

    # degrees-of-freedom (dof) per node
    ndof = 1

    # Number of element nodes (bilinear elements)
    nen = 4

    # Number of equations
    neq = nnp*ndof

    return nel, lpx, lpy, nnp, ndof, nen, neq


def phys_coord(nelx, nely):
    """ This function returns the physical coordinates of the nodes.

    Input:
    ------
    nelx:   integer
            number of elements in the x direction.
    nely:   integer
            number of elements in the y direction.

    Output:
    -------
    x:      float (1d array)
            the coordinate of the node in the x direction
    y:      float (1d array)
            the coordinate of the node in the y direction

    The geometry we are working on is like the following.
    (for nelx = 2, nely = 2)
    6---------7----------8
    |         |   (3)    |
    |   (2)   |      ----5
    |      ---4-----/    |
    3-----/   |   (1)    |
    |         |      ----2
    |   (0)   |     /
    |     ----1----/
    0----/
    There are 4 elements (numbering in parenthesis), and 9 nodes.
    Bottom edge (0 to 2) is y=0.5x^2. (see src/test_subroutines.py)
    This function returns x,y as 9x2 array for the above mesh.
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Divide [0,1] by lpx (mesh in the x direction)
    x0 = np.linspace(0, 1, lpx)
    y0 = 0.5 * x0**2               # the bottom geometry line

    y = np.zeros((nnp, 1))
    for i in range(0, lpx):
        # Divide [0,1] by lpy (mesh in the y direction)
        y1 = np.linspace(y0[i], 1, lpy)
        for j in range(0, lpy):
            y[i + j*lpx] = y1[j]   # collection of y

    x = np.zeros((nnp, 1))
    for i in range(0, lpy):
        for j in range(0, lpx):
            x[j + i*lpx] = x0[j]   # collection of x

    return x, y


def connectivity(nelx, nely):
    """ This function returns the connectivity matrix.

    Input:
    ------
    nelx:   integer
            number of elements in the x direction
    nely:   integer
            number of elements in the y direction

    Output:
    ------
    A:      integer (2d array)
            connectivity matrix, IEN
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Total number of elements in the domain
    nel = nelx*nely
    IEN = np.zeros((nen, nel), dtype=int)
    rowcount = 0
    # Connectivity matrix for 4-node elements
    if nen == 4:
        for elementcount in range(0, nel):
            IEN[0][elementcount] = elementcount + rowcount
            IEN[1][elementcount] = elementcount + 1 + rowcount
            IEN[2][elementcount] = elementcount + (lpx + 1) + rowcount
            IEN[3][elementcount] = elementcount + (lpx) + rowcount
            if np.mod(elementcount + 1, lpx - 1) == 0:
                rowcount = rowcount + 1

    return IEN


def Dirichlet_BCs(nelx, nely, T0_bottom, T0_left):
    """ This function returns the flags and e_bc arrays.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    T0_bottom:  float
                prescribed temperature on the bottom edge
    T0_left:    float
                prescribed temperature on the left edge

    Output:
    -------
    e_bc:       float (1d array)
                contains the prescribed temperature on the boundary
    flags:      boolean (1d array)
                indicates nodes that have prescribed temperature
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Array to set B.C. flags
    flags = np.zeros((neq, 1), dtype=int)

    # Essential B.C. array
    e_bc = np.zeros((neq, 1))

    # Essential B.C. (prescribed temperature)
    for i in range(0, lpx):
        flags[i] = 1
        e_bc[i] = T0_bottom       # bottom edge

    for i in range(lpx, nnp - nelx, lpx):
        flags[i] = 1
        e_bc[i] = T0_left         # left edge

    return e_bc, flags


def setup_ID_LM(nelx, nely, T0_bottom, T0_left):
    """ This function returns the ID array and LM matrix.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    T0_bottom:  float
                prescribed temperature on the bottom edge
    T0_left:    float
                prescribed temperature on the left edge

    Output:
    -------
    d0:         integer (1d array)
                dof for each node
    ID:         integer (2d array)
                mapping between each element and their nodes
    LM:         integer
                dof for each element of ID
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Number of nodes on essential boundaries (left and bottom edges)
    nd = lpx + lpy - 1

    # Array of node numbers
    ID = np.zeros((neq, 1), dtype=int)

    # Initialize nodal temperature vector
    d0 = np.zeros((neq, 1))

    # Read connectivity matrix, IEN
    IEN = connectivity(nelx, nely)

    # Read prescribed temperature on the boundaries and flags
    e_bc, flags = Dirichlet_BCs(nelx, nely, T0_bottom, T0_left)
    LM = np.zeros((4, nel), dtype=int)
    count = 0
    count1 = 0
    for i in range(0, neq):
        if flags[i] == 1:
            ID[i][0] = count
            d0[count] = e_bc[i]
            count = count + 1
        else:
            ID[i][0] = count1 + nd
            count1 = count1 + 1

    for i in range(0, nel):
        for j in range(0, nen):
            LM[j][i] = ID[IEN[j][i]][0]

    return d0, ID, LM


def basis(xi, eta):
    """ This function returns shape functions for a 4 node bilinear element
    After reading the coordinate of each element in (x,y), (x,y) is mapped
    to (xi, eta) for the integration.

    Input:
    ------
    xi, eta:    float
                natural coordinate system in the range [-1, 1]

    Output:
    ------
    N:          float (1d array)
                N = [N0, N1, N2, N3]
                defined on the natural coordinates (xi, eta)
    """
    N = 0.25 * np.array([[(1-xi)*(1-eta),
                          (1+xi)*(1-eta),
                          (1+xi)*(1+eta),
                          (1-xi)*(1+eta)]])
    return N


def d_basis(xi, eta, coord):
    """ This function returns the derivative of shape functions
    for a 4 node bilinear element. (x, y) is used to compute
    Jacobian for maping to (xi, eta).

    Input:
    ------
    xi, eta:    float
                natural coordinate system in the range [-1, 1]
    coord:      float
                physical coordinates of elements in (x, y)

    Output:
    ------
    B:          float (2d array)
                derivative of shape functions
                B=[dN/dx; dN/dy]
    detJ:       float (2x4 array)
                determinant of B (the Jacobian matrix)

    """
    # Derivative of N
    #   the first row with respect to xi
    #   the second row with respect to eta
    dN = 0.25*np.array([[eta-1, 1-eta, 1+eta, -eta-1],
                        [xi-1, -xi-1, 1+xi, 1-xi]])

    # Compute Jacobian matrix
    J = np.matmul(dN, coord)

    # Determinate of the Jacobian matrix
    detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    invJ = np.zeros((2, 2))
    invJ[0][0] = J[1][1]
    invJ[0][1] = -J[0][1]
    invJ[1][0] = -J[1][0]
    invJ[1][1] = J[0][0]

    # Inverse of the Jacobian matrix
    invJ = invJ/detJ

    # Derivative of N in (x,y) coordinates
    #   B[0, :] = dN/dx
    #   B[1, :] = dN/dy
    B = np.matmul(invJ, dN)

    return B, detJ


def gauss(ngp):
    """ This function returns the quadrature weight and the quadrature points.

    Input:
    ------
    ngp:        integer
                number of Gauss points for the integration
                available choices are only 1 and 2

    Output:
    ------
    w:          float
                quadrature weight for each quadrature point
    gp:         float
                quadrature coordinates in the natural coordinate
                system (xi,eta) for each quadrature point
    """
    if ngp == 1:
        gp = 0
        w = 2
        return w, gp
    elif ngp == 2:
        gp = np.array([-0.57735027, 0.57735027])
        w = np.array([1, 1])
        return w, gp
    else:
        print("Error: This code supports only 2 quadrature points so far!")
        sys.exit(1)

    return None


def heat2delem(nelx, nely, ngp, s0, e):
    """ This function returns the stiffness, ke, and forcing function, fe.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    ngp:        integer
                number of Gauss points for the integration
    e:          integer
                element number
    s0:         float
                magnitude of the heat source

    Output:
    ------
    ke:         float (4x4 matrix)
                stiffness matrix
    fe:         float (4x1 vector)
                forcing vector
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Initialize element conductivity matrix
    ke = np.zeros((nen, nen))

    # Initialize element nodal source vector
    fe = np.zeros((nen, 1))

    # Get coordinates of element nodes
    je = np.zeros((nel, 1), dtype=int)
    IEN = connectivity(nelx, nely)
    for i in range(nen):
        # For each element e, we save the node numbers given in IEN
        je[i][0] = IEN[i][e]

    # Read the (x,y) coordinates of all elements
    x, y = phys_coord(nelx, nely)

    # Find the (x,y) coordinates of element e based on je above
    C = np.array([[x[je[0][0]][0],
                   x[je[1][0]][0],
                   x[je[2][0]][0],
                   x[je[3][0]][0]],
                  [y[je[0][0]][0],
                   y[je[1][0]][0],
                   y[je[2][0]][0],
                   y[je[3][0]][0]]])
    C = np.transpose(C)

    # Get gauss points and weights
    w, gp = gauss(ngp)

    # Material properties
    k = 5                        # Thermal conductivity
    D = k*np.identity(2)         # Conductivity matrix
    s = s0*np.ones((nen, nel))   # Heat source

    # Compute element conductivity matrix and nodal flux vector
    # Loop over Gauss points to compute the integral
    for i in range(ngp):
        for j in range(ngp):
            # Get reference coordinates
            eta = gp[i]
            xi = gp[j]

            # Shape functions matrix
            N = basis(xi, eta)

            # Derivative of the shape functions
            B, detJ = d_basis(xi, eta, C)

            # Element conductivity matrix
            wwdetJ = w[i] * w[j] * detJ
            ke = ke + wwdetJ * np.matmul(np.matmul(np.transpose(B), D), B)

            # Compute s(x)
            S = np.array([[s[0][e]], [s[1][e]],
                          [s[2][e]], [s[3][e]]])
            se = np.matmul(N, S)

            # Element nodal source vector
            fe = fe + w[i] * w[j] * np.matmul(np.transpose(N), se) * detJ

    return ke, fe


def assembly(nelx, nely, ngp, s0, T0_bottom, T0_left):
    """ This function assembles stiffness, K, and forcing, F.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    ngp:        integer
                number of Gauss points for the integration
    s0:         float
                magnitude of the heat source
    T0_bottom:  float
                prescribed temperature on the bottom edge
    T0_left:    float
                prescribed temperature on the left edge

    Output:
    ------
    K:          float (2d array)
                assembled stiffness matrix
    F:          float (1d array)
                assembled forcing vector
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Get LM
    d0, ID, LM = setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    f = np.zeros((neq, 1))         # Initialize global forcing vector
    K = np.zeros((neq, neq))       # Initialize global stiffness matrix
    for e in range(nel):
        # Get element stiffness and force for each element
        ke, fe = heat2delem(nelx, nely, ngp, s0, e)
        for loop1 in range(nen):
            i = LM[loop1][e]
            # Assemble the forcing vector
            f[i] = f[i] + fe[loop1]
            for loop2 in range(nen):
                j = LM[loop2][e]
                # Assemble the stiffness matrix
                K[i][j] = K[i][j] + ke[loop1][loop2]

    return K, f


def NeumannBCs(nelx, nely, flux_top):
    """ This function returns the number of elements
    on top boundary with Neumann BCs and array n_bc
    which contains the prescribed flux on top edge

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    flux_top:   float
                prescribed flux on the top edge

    Output:
    -------
    nbe:        integer
                number of elements with non-zero flux on the top boundary
    n_bc:       float()
                n_bc[0,:] => node number
                n_bc[0,:] => prescribed flux for each node on the top boundary
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # natural B.C. on the top edge
    n_bc = np.zeros((2, lpx))
    for i in range(0, lpx):
        # Nodes' numbers on the top edge
        n_bc[0][i] = nnp-lpx + i
        # Flux values on the corresponding nodes
        n_bc[1][i] = flux_top

    # Number of elements with non-zero flux on the boundary
    nbe = nelx

    return nbe, n_bc


def src_flux(nelx, nely, T0_bottom, T0_left, flux_top, ngp, F):
    """ This function computes the flux on the top boundary
    and adds it to the global forcing vector, F.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    T0_bottom:  float
                prescribed temperature on the bottom edge
    T0_left:    float
                prescribed temperature on the left edge
    flux_top:   float
                prescribed flux on the top edge
    ngp:        integer
                number of Gauss points for the integration
    F:          float (1d array)
                global forcing vector

    Output:
    ------
    F:          float (1d array)
                global forcing vector (updated)
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Initialize point source defined at a node
    P = np.zeros((neq, 1))
    d0, ID, LM = setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    nbe, n_bc = NeumannBCs(nelx, nely, flux_top)

    # Updated F, if there is an external source on the nodes
    for i in range(0, neq):
        F[ID[i]] = F[ID[i]] + P[ID[i]]

    x, y = phys_coord(nelx, nely)
    for i in range(0, nbe):
        fq = np.zeros((2, 1))
        n_bce = np.zeros((2, 1))

        node1 = int(n_bc[0][i])
        node2 = int(n_bc[0][i+1])
        n_bce[0][0] = n_bc[1][i]
        n_bce[1][0] = n_bc[1][i+1]

        x1 = x[node1]
        y1 = y[node1]
        x2 = x[node2]
        y2 = y[node2]

        length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
        J = length/2

        w, gp = gauss(ngp)

        for i in range(0, ngp):
            N = np.zeros((2, 1))
            psi = gp[i]
            N[0][0] = 0.5*(1-psi)
            N[1][0] = 0.5*(1+psi)

            flux = N[0][0]*n_bce[0][0] + N[1][0]*n_bce[1][0]
            fq[0][0] = fq[0][0] + w[i]*N[0][0]*flux*J
            fq[1][0] = fq[1][0] + w[i]*N[1][0]*flux*J

        fq = -fq
        # Update F by adding flux on the top boundary
        F[ID[node1]] = F[ID[node1]] + fq[0][0]
        F[ID[node2]] = F[ID[node2]] + fq[1][0]

    return F


def solvedr(nelx, nely, K, F, d0):
    """ This function solves K*d=F and returns d.

    Input:
    ------
    nelx:       integer
                number of elements in the x direction
    nely:       integer
                number of elements in the y direction
    K:          float (2d array)
                global stiffness matrix
    F:          float (1d array)
                global forcing vector
    d0:         float
                initialzed temperature

    Output:
    ------
    d:          float (1d array)
                the final solution
                temperature at all nodes in the domain
    """
    # Get the setup properties
    nel, lpx, lpy, nnp, ndof, nen, neq = setup(nelx, nely)

    # Number of nodes on the essential boundaries (left and bottom edges)
    nd = lpx+lpy - 1

    K_E = np.zeros((nd, nd))
    for i in range(nd):
        for j in range(nd):
            K_E[i][j] = K[i][j]

    K_F = np.zeros((neq-nd, neq-nd))
    for i in range(neq-nd):
        for j in range(neq-nd):
            K_F[i][j] = K[nd+i][nd+j]

    K_EF = np.zeros((nd, neq-nd))
    for i in range(nd):
        for j in range(neq-nd):
            K_EF[i][j] = K[i][nd+j]
    K_EF = np.transpose(K_EF)

    f_E = np.zeros((neq-nd, 1))
    for i in range(neq-nd):
        f_E[i][0] = F[nd+i][0]

    d_E = np.zeros((nd, 1))
    for i in range(nd):
        d_E[i][0] = d0[i][0]

    d_F = np.zeros((neq-nd, 1))
    d_F = np.linalg.solve(K_F, f_E - np.matmul(K_EF, d_E))

    d = np.zeros((neq, 1))
    for i in range(nd):
        d[i][0] = d_E[i][0]

    for i in range(neq-nd):
        d[i+nd][0] = d_F[i][0]

    return d
