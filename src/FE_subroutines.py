import numpy as np
import math


def phys_coord(nelx, nely):
    """ This function return the physical coordinates of the nodes.
    Input:
    ------
    nelx: number of elements in x direction.
    nely: number of elements in y direction.

    Return:
    -------
    x,y the coordinates of each nodes.
    The geometry we are working on is like
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
    We have 4 elements (numbering in parenthesis ), and 9 nodes.
    Bottom edge (0 to 2) is y=0.5x^2. (see unit test)
    So this functions return x,y as 9x2 array for above mesh.
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # divide [0,1] by lpx (mesh in x direction)
    x0 = np.linspace(0, 1, lpx)
    y0 = 0.5 * x0**2               # the bottom geometry line

    y = np.zeros((nnp, 1))
    for i in range(0, lpx):
        # divide [0,1] by lpy (mesh in y direction)
        y1 = np.linspace(y0[i], 1, lpy)
        for j in range(0, lpy):
            y[i + j*lpx] = y1[j]   # collection of y coordinate

    x = np.zeros((nnp, 1))
    for i in range(0, lpy):
        for j in range(0, lpx):
            x[j + i*lpx] = x0[j]   # collection of x coordinate
    return x, y


def connectivity(nelx, nely):
    """ This function returns the connectivity matrix.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction

    Return:
    ------
    A connectivity matrix, IEN
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # total number of elements in a domain
    nel = nelx*nely
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


def DirichletBCs(nelx, nely, T0_bottom, T0_left):
    """ This function return the flags and e_bc arrays.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    T0_bottom: prescribed temperature on the bottom edge
    T0_left: prescribed temperature on the left edge

    Return:
    -------
    e_bc: array that contains the prescribed temperature on the boundary
    flags: Boolean array to find which nodes has prescribed temperature
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # degrees-of-freedom per node
    ndof = 1
    # number of equation
    neq = nnp*ndof
    # array to set B.C flags
    flags = np.zeros((neq, 1), dtype=int)
    # essential B.C array
    e_bc = np.zeros((neq, 1))
    # essential B.C. (prescribed temperature)
    for i in range(0, lpx):
        flags[i] = 2
        e_bc[i] = T0_bottom        # bottom edge

    for i in range(lpx, nnp - nelx, lpx):
        flags[i] = 2
        e_bc[i] = T0_left         # left edges

    return e_bc, flags


def setup_ID_LM(nelx, nely, T0_bottom, T0_left):
    """ This function return the ID array and LM matrix.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    T0_bottom: prescribed temperature on the bottom edge
    T0_left: prescribed temperature on the left edge

    Return:
    -------
    ID and LM matrices
    """
    # total number of elements in a domain
    nel = nelx*nely
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # degrees-of-freedom per node
    ndof = 1
    # number of element nodes (element is bilinear)
    nen = 4
    # number of equation
    neq = nnp*ndof
    # number of nodes on essential boundary (left+bottom edge)
    nd = lpx+lpy - 1
    # array the contains node numbers
    ID = np.zeros((neq, 1), dtype=int)
    # initialize nodal temperature vector
    d0 = np.zeros((neq, 1))
    # read connectivity matrix, IEN
    IEN = connectivity(nelx, nely)
    # read, prescribed temperature on the boundary and flags
    e_bc, flags = DirichletBCs(nelx, nely, T0_bottom, T0_left)
    LM = np.zeros((4, nel), dtype=int)
    count = 0
    count1 = 0
    for i in range(0, neq):
        if flags[i] == 2:
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
    """ This function returns the shape function for a 4 nodes bilinear element

    Input:
    ------
    xi, eta are the natural coordinate system.
    natural coordinate system is always in the range [-1, 1]
    So, after reading the coordinate of each element in (x,y)
    we map (x,y) to (xi, eta) for integration.

    Return:
    ------
    an array of shape function as N = [N1, N2, N3, N4]
    which defined on the natural coordinate (xi,eta)
    """
    N = 0.25 * np.array([[(1-xi)*(1-eta),
                          (1+xi)*(1-eta),
                          (1+xi)*(1+eta),
                          (1-xi)*(1+eta)]])
    return N


def d_basis(xi, eta, coord):
    """ This function returns the  derivative of shape function
    for a 4 nodes bilinear element.

    Input:
    ------
    xi, eta are the natural coordinate system.
    coord: physical coordinate of element in (x, y)
    we use (x, y) to compute Jacobian for maping to (xi, eta)

    Return:
    ------
    Derivative of shape function B=[dN/dx; dN/dy]
    """
    # derivative of N, first row respect to xi,
    # and second row of N respect to eta
    dN = 0.25*np.array([[eta-1, 1-eta, 1+eta, -eta-1],
                        [xi-1, -xi-1, 1+xi, 1-xi]])
    # compute Jacobian matrix
    J = np.matmul(dN, coord)
    # determinate of Jacobian
    detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    invJ = np.zeros((2, 2))
    invJ[0][0] = J[1][1]
    invJ[0][1] = -J[0][1]
    invJ[1][0] = -J[1][0]
    invJ[1][1] = J[0][0]
    # Inverse of Jacobian
    invJ = invJ/detJ
    # derivative of N in (x,y) coordinate
    # first row of B is dN/dx
    # second row of B is dN/dy
    # size B is 2x4
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
    gp: quadrature coordinates in the natural coordinate system (xi,eta)
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


def heat2delem(nelx, nely, ngp, s0, e):
    """ This function returns the stiffness, ke, and forcing function, fe.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    ngp: number of Gauss point you want for integration
    e: element number
    s0: magnitude of the heat source

    Return:
    ------
    ke: a 4x4 stiffness matrix
    fe: a 4x1 forcing vector
    """
    # total number of elements in a domain
    nel = nelx*nely
    # number of element nodes
    nen = 4
    # Initialize element conductance matrix
    ke = np.zeros((nen, nen))
    # Initialize element nodal source vector
    fe = np.zeros((nen, 1))

    # Get coordinates of element nodes
    je = np.zeros((nel, 1), dtype=int)
    IEN = connectivity(nelx, nely)
    for i in range(nen):
        # for each element e, we save the node numbers given in IEN
        je[i][0] = IEN[i][e]
    # read the (x,y) coordinates of all elements
    x, y = physCoord(nelx, nely)
    # find the (x,y) coordinate of element e based on je above
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
    # material property
    k = 5                        # thermal conductivity
    D = k*np.identity(2)         # conductivity matrix
    s = s0*np.ones((nen, nel))   # heat source
    # compute element conductance matrix and nodal flux vector
    # this part is loop over Gauss point to compute integral
    for i in range(ngp):
        for j in range(ngp):
            # Get reference coordinates
            eta = gp[i]
            xi = gp[j]

            # Shape functions matrix
            N = basis(xi, eta)

            # Derivative of the shape functions
            B, detJ = d_basis(xi, eta, C)

            # element conductance matrix
            wwdetJ = w[i] * w[j] * detJ
            ke = ke + wwdetJ * np.matmul(np.matmul(np.transpose(B), D), B)

            # compute s(x)
            S = np.array([[s[0][e]], [s[1][e]],
                          [s[2][e]], [s[3][e]]])
            se = np.matmul(N, S)

            # element nodal source vector
            fe = fe + w[i] * w[j] * np.matmul(np.transpose(N), se) * detJ

    return ke, fe


def assembly(nelx, nely, ngp, s0, T0_bottom, T0_left):
    """ This function assembles stiffness, K, and forcing, F.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    ngp: number of Gauss point you want for integration
    s0: magnitude of the heat source
    T0_bottom: prescribed temperature on the bottom edge
    T0_left: prescribed temperature on the left edge

    Return:
    ------
    K: assembled stiffness matrix
    F: assembled forcing vector
    """
    # total number of elements in a domain
    nel = nelx*nely
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # degrees-of-freedom per node
    ndof = 1
    # number of element nodes (element is bilinear)
    nen = 4
    # number of equation
    neq = nnp*ndof
    # get LM arrays
    d0, ID, LM = setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    f = np.zeros((neq, 1))         # initialize global force vector
    K = np.zeros((neq, neq))       # initialize global stiffness matrix
    for e in range(nel):
        # Get element stiffness and force for each element
        ke, fe = heat2delem(nelx, nely, ngp, s0, e)
        for loop1 in range(nen):
            i = LM[loop1][e]
            f[i] = f[i] + fe[loop1]  # Assemble forces
            for loop2 in range(nen):
                j = LM[loop2][e]
                K[i][j] = K[i][j] + ke[loop1][loop2]  # Assemble stiffness

    return K, f


def NeumannBCs(nelx, nely, flux_top):
    """ This function return the number of elements
    on top boundary with Neumann BCs and array n_bc
    which contains the prescribed flux on top edge

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    flux_top: prescribed flux on the top edge

    Return:
    -------
    n_bc: array that contains the prescribed flux on the top boundary
    nbe: number of element with non-zero flux on the top boundary
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # natural B.C  - on top edge,
    n_bc = np.zeros((2, lpx))
    for i in range(0, lpx):
        # node's numbers on top
        n_bc[0][i] = nnp-lpx + i
        # flux on those nodes
        n_bc[1][i] = flux_top
    # number of element with non-zero flux on the boundary
    nbe = nelx

    return nbe, n_bc


def src_flux(nelx, nely, T0_bottom, T0_left, flux_top, ngp, F):
    """ This function computes the flux on the top boundary
    then adds to the global forcing vector, F, we have
    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    T0_bottom: prescribed temperature on the bottom edge
    T0_left: prescribed temperature on the left edge
    flux_top: prescribed flux on the top edge
    ngp: number of Gauss point you want for integration
    F: global forcing vector

    Return:
    ------
    updated global forcing vector
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # degrees-of-freedom per node
    ndof = 1
    # number of equation
    neq = nnp*ndof
    # initialize point source defined at a node
    P = np.zeros((neq, 1))
    d0, ID, LM = setup_ID_LM(nelx, nely, T0_bottom, T0_left)
    nbe, n_bc = NeumannBCs(nelx, nely, flux_top)
    # updated F, if we have external source on the nodes
    for i in range(0, neq):
        F[ID[i]] = F[ID[i]] + P[ID[i]]

    x, y = physCoord(nelx, nely)
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
        # update F by adding flux on the top boundary
        F[ID[node1]] = F[ID[node1]] + fq[0][0]
        F[ID[node2]] = F[ID[node2]] + fq[1][0]
    return F


def solvedr(nelx, nely, K, F, d0):
    """ This function solves K*d=F and returns the d.

    Input:
    ------
    nelx: number of elements in x direction
    nely: number of elements in y direction
    K: global stiffness matrix
    F: global forcig vector
    d0: initialzed temperature

    Return:
    ------
    d: an array of the solution which is the temperatures
    at all nodes we have in the domain
    """
    # number of nodes in x direction
    lpx = nelx + 1
    # number of nodes in y direction
    lpy = nely + 1
    # total number of nodes in a domain
    nnp = lpx*lpy
    # degrees-of-freedom per node
    ndof = 1
    # number of equation
    neq = nnp*ndof
    # number of nodes on essential boundary (left+bottom edge)
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
