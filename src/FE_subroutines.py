import numpy as np
import variables as var
import math


def physCoord(lpx, lpy):
    """ This function return the physical coordinates of the nodes.
    Input:
    ------
    lpx: number of nodes in x direction.
    lpy: number of nodes in y direction.

    Return:
    -------
    x,y the coordinates of each nodes
    """

    x0 = np.linspace(0, 1, lpx)
    y0 = 0.5 * x0**2               # the bottom geometry line

    y = np.zeros((var.nnp, 1))
    for i in range(0, lpx):
        y1 = np.linspace(y0[i], 1, lpy)
        for j in range(0, lpy):
            y[i + j*lpx] = y1[j]   # collection of y coordinate

    x = np.zeros((var.nnp, 1))
    for i in range(0, lpy):
        for j in range(0, lpx):
            x[j + i*lpx] = x0[j]   # collection of x coordinate
    return x, y


# generate the IEN connectivity array
def connectivity(nel, lpx):
    """ This function returns the connectivity matrix.

    Input:
    ------
    nel: number of elements in a domain
    lpx: number of nodes in x direction

    Return:
    ------
    A connectivity matrix, IEN
    """

    IEN = np.zeros((4, nel), dtype=int)
    rowcount = 0
    for elementcount in range(0, nel):
        IEN[0][elementcount] = elementcount + rowcount
        IEN[1][elementcount] = elementcount + 1 + rowcount
        IEN[2][elementcount] = elementcount + (lpx + 1) + rowcount
        IEN[3][elementcount] = elementcount + (lpx) + rowcount
        if np.mod(elementcount + 1, lpx - 1) == 0:
            rowcount = rowcount + 1

    return IEN


def setup_ID_LM(neq, nel, nd):
    """ This function return the ID array and LM matrix.

    Input:
    ------
    neq: the number of equations we will solve
    nel: the number of elements in the domain
    nd: number of nodes on the boundary with Dirichlet condition

    Return:
    -------
    ID and LM matrices
    """

    IEN = connectivity(nel, var.lpx)
    LM = np.zeros((4, nel), dtype=int)
    count = 0
    count1 = 0
    for i in range(0, neq):
        if var.flags[i] == 2:
            var.ID[i] = count
            var.d[count] = var.e_bc[i]
            count = count + 1
        else:
            var.ID[i] = count1 + nd
            count1 = count1 + 1

    for i in range(0, nel):
        for j in range(0, var.nen):
            LM[j][i] = var.ID[IEN[j][i]]

    return LM


def basis(xi, eta):
    """ This function returns the shape function for a 4 nodes bilinear element

    Input:
    ------
    xi, eta are the natural coordinate system

    Return:
    ------
    an array of shape function as N = [N1, N2, N3, N4]
    """
    N = 0.25 * np.array([[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]])
    return N


def d_basis(xi, eta, coord):
    """ This function returns the  derivative of shape function
    for a 4 nodes bilinear element.

    Input:
    ------
    xi, eta are the natural coordinate system.
    coord: coordinates of the node to compute Jacobian

    Return:
    ------
    Derivative of shape function B=[dN/dx; dN/dy]
    """
    # Calculate the Grad(N) matrix
    dN = 0.25*np.array([[eta-1, 1-eta, 1+eta, -eta-1],
                        [xi-1, -xi-1, 1+xi, 1-xi]])

    J = np.matmul(dN, coord)      # compute Jacobian matrix

    detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    invJ = np.zeros((2, 2))
    invJ[0][0] = J[1][1]
    invJ[0][1] = -J[0][1]
    invJ[1][0] = -J[1][0]
    invJ[1][1] = J[0][0]
    invJ = invJ/detJ
    B = np.matmul(invJ, dN)
    return B, detJ


def gauss(ngp):
    """ This function returns quadrature weight and quadrature points

    Input:
    ------
    ngp: number of Gauss point you want for integration

    Return:
    ------
    w: weight of quadrature
    gp: quadrature coordinates in the natural coordinate system (xi,eta) => in
    [-1, 1] domain
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
        print("Error: This code supports only 1 and 2 quadrature points.")


def heat2delem(e):
    """ This function returns the stiffness, ke, and forcing function, fe.

    Input:
    ------
    e: element number

    Return:
    ------
    ke: a 4x4 stiffness matrix
    fe: a 4x1 forcing vector
    """
    ke = np.zeros((var.nen, var.nen))  # Initialize element conductance matrix
    fe = np.zeros((var.nen, 1))      # Initialize element nodal source vector

    # Get coordinates of element nodes
    je = np.zeros((var.nel, 1), dtype=int)
    IEN = connectivity(var.nel, var.lpx)
    for i in range(var.nen):
        je[i][0] = IEN[i][e]

    x, y = physCoord(var.lpx, var.lpy)
    C = np.array([[x[je[0][0]][0], x[je[1][0]][0], x[je[2][0]][0], x[je[3][0]][0]],
                 [y[je[0][0]][0], y[je[1][0]][0], y[je[2][0]][0], y[je[3][0]][0]]])
    C = np.transpose(C)

    # Get gauss points and weights
    w, gp = gauss(var.ngp)

    # compute element conductance matrix and nodal flux vector 
    for i in range(var.ngp):
        for j in range(var.ngp):
            # Get reference coordinates
            eta = gp[i]           
            xi = gp[j]

            # Shape functions matrix
            N = basis(xi, eta)
            # Derivative of the shape functions
            B, detJ = d_basis(xi, eta, C)
            # element conductance matrix
            ke = ke + w[i] * w[j] * np.matmul(np.matmul(np.transpose(B), var.D), B) * detJ

            # compute s(x)
            s = np.array([[var.s[0][e]], [var.s[1][e]],
            [var.s[2][e]], [var.s[3][e]]])
            se = np.matmul(N, s)

            # element nodal source vector
            fe = fe + w[i] * w[j] * np.matmul(np.transpose(N), se) * detJ

    return ke, fe


def assembly(e):
    """ This function the assembled stiffness, K, and forcing, F.

    Input:
    ------
    e: element number

    Return:
    ------
    K: assembled stiffness matrix
    F: assembled forcing vector
    """
    # Get element stiffness and force
    ke, fe = heat2delem(e)
    LM = setup_ID_LM(var.neq, var.nel, var.nd)
    for loop1 in range(var.nen):
        i = LM[loop1][e]
        var.f[i] = var.f[i] + fe[loop1]  # Assemble forces
        for loop2 in range(var.nen):
            j = LM[loop2][e]
            var.K[i][j] = var.K[i][j] + ke[loop1][loop2]  # Assemble stiffness
    return


def src_flux(neq, nbe, ngp):
    """ This function computes the flux on the boundary
    that we have non-zero flux
    Input:
    ------
    neq: number of equations
    nbe: number of element on the boundary with non-zero flux
    ngp: number of Gauss point for integration

    Return:
    ------
    a forcing flux vector
    """
    for i in range(0, neq):
        var.f[var.ID[i]] = var.f[var.ID[i]] + var.P[var.ID[i]]

    x, y = physCoord(var.lpx, var.lpy)
    for i in range(0, nbe):
        fq = np.zeros((2, 1))
        n_bce = np.zeros((2, 1))

        node1 = int(var.n_bc[0][i])
        node2 = int(var.n_bc[0][i+1])
        n_bce[0][0] = var.n_bc[1][i]
        n_bce[1][0] = var.n_bc[1][i+1]

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

        var.f[var.ID[node1]] = var.f[var.ID[node1]] + fq[0][0]
        var.f[var.ID[node2]] = var.f[var.ID[node2]] + fq[1][0]
    return


def solvedr(neq, nd):
    """ This function solves K*d=F and returns the d.

    Input:
    ------
    neq: number of equations
    nd: number of nodes on the boundary with Dirichlet condition

    Return:
    ------
    d: an array of the solution which is the temperatures
    at all nodes we have in the domain
    """
    K_E = np.zeros((nd, nd))
    for i in range(nd):
        for j in range(nd):
            K_E[i][j] = var.K[i][j]

    K_F = np.zeros((neq-nd, neq-nd))
    for i in range(neq-nd):
        for j in range(neq-nd):
            K_F[i][j] = var.K[nd+i][nd+j]

    K_EF = np.zeros((nd, neq-nd))
    for i in range(nd):
        for j in range(neq-nd):
            K_EF[i][j] = var.K[i][nd+j]
    K_EF = np.transpose(K_EF)

    f_E = np.zeros((neq-nd, 1))
    for i in range(neq-nd):
        f_E[i][0] = var.f[nd+i][0]

    d_E = np.zeros((nd, 1))
    for i in range(nd):
        d_E[i][0] = var.d[i][0]

    d_F = np.zeros((neq-nd, 1))
    d_F = np.linalg.solve(K_F, f_E - np.matmul(K_EF, d_E))

    d = np.zeros((neq, 1))
    for i in range(nd):
        d[i][0] = d_E[i][0]

    for i in range(neq-nd):
        d[i+nd][0] = d_F[i][0]

    return d
