import numpy as np
import variables as var
import math

#=============================================================
#               Mesh Generation and IEN (mesh2d.m)
#=============================================================

# This function returns the physical coordinates of the nodes
def physCoord(lpx, lpy):

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

    IEN = np.zeros((4, nel), dtype = int)
    rowcount = 0
    for elementcount in range(0,nel):
        IEN[0][elementcount] = elementcount + rowcount
        IEN[1][elementcount] = elementcount + 1 + rowcount
        IEN[2][elementcount] = elementcount + (lpx + 1) + rowcount
        IEN[3][elementcount] = elementcount + (lpx) + rowcount
        if np.mod(elementcount+1,lpx-1) == 0:
            rowcount = rowcount + 1
            
    return IEN

#=============================================================
#               Plot Mesh (plotmesh.m)
#=============================================================

 


#=============================================================
#               Setup ID and LM (setup_ID_LM.m)
#=============================================================
def setup_ID_LM (neq, nel, nd):

    IEN = connectivity(nel, var.lpx)
    LM = np.zeros((4, nel), dtype = int)
    count  = 0
    count1 = 0
    for i in range(0, neq):
        if var.flags[i] == 2:
            var.ID[i]  = count
            var.d[count] = var.e_bc[i]
            count  = count +1
        else:
            var.ID[i]  = count1 + nd
            count1 = count1 + 1
  
    for i in range(0, nel):
        for j in range(0, var.nen):
            LM[j][i] = var.ID[IEN[j][i]]

    return LM

#=============================================================
#               basis and derivative of basis and gauss
#=============================================================

# N defines bilinear basis function
# N = [N1, N2, N3, N4] 
def basis(xi, eta):
    N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
    return N


def d_basis(xi, eta, coord):

    #Calculate the Grad(N) matrix
    dN = np.array(0.25*[[eta-1, 1-eta, 1+eta,-eta-1],
                        [xi-1 , -xi-1, 1+xi , 1-xi]])

    J     = dN*coord      # compute Jacobian matrix 
    #detJ  = np.linalg.det(J)     # Jacobian  
    #B     = np.linalg.solve(J, dN)       # compute the B matrix

    detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0]
    invJ = np.zeros((2,2))
    invJ[0][0] =  J[1][1]
    invJ[0][1] = -J[0][1]
    invJ[1][0] = -J[1][0]
    invJ[1][1] =  J[0][0]
    invJ = invJ/detJ
    B     = invJ*dN
    return B, detJ


# nqpts: number of quadrature points
# qpts:  coordinates of the quadrature points
# w:     quadrature weight
def gauss(ngp):

    if ngp == 1:
        gp = 0
        w = 2
        return gp, w
    elif ngp == 2:
        gp = [-0.57735027,  0.57735027 ]  
        w  = [1,            1          ]
        return gp, w
    else:
        print("Error: This code supports only 1 and 2 quadrature points.")

    return
    


#=============================================================
#               stiffness, force element (heat2Delem.m)
#=============================================================
# e: number of elements
def heat2delem(e):

    ke = np.zeros(var.nen, var.nen) # Initialize element conductance matrix
    fe = np.zeros(var.nen, 1)      # Initialize element nodal source vector

    # Get coordinates of element nodes
    je = np.zeros((var.nel, 1))
    IEN = connectivity(var.nel, var.lpx)
    for i in range(var.nen):
        je[i] = IEN[i][e]

    x, y = physCoord(var.lpx, var.lpy)
    # TODO: make the size of x,y dynamic
    C  = np.transpose( [ [ x[je[0]], x[je[1]], x[je[2]], x[je[3]] ] , 
                         [ y[je[0]], y[je[1]], y[je[2]], y[je[3]] ] 
                       ] )

    # Get gauss points and weights
    w, gp = gauss(var.ngp)

    # compute element conductance matrix and nodal flux vector 
    for i in range(var.ngp):
        for j in range(var.ngp):
            # Get reference coordinates
            eta = gp[i]           
            psi = gp[j]

            # Shape functions matrix
            N = basis(eta, psi)

            # Derivative of the shape functions 
            B, detJ = d_basis(eta, psi, C)

            # element conductance matrix
            ke = ke + w[i] * w[j] * np.transpose(B) * var.D * B * detJ

            # compute s(x)
            s = [ var.s[0][e], var.s[1][e], var.s[2][e], var.s[3][e] ]
            se = N * s

            # element nodal source vector
            fe = fe + w[i] * w[j] * np.transpose(N) * se * detJ

    return ke, fe

#=============================================================
#               Assembly (assembly.m)
#=============================================================
# Assemble element matrices and vectors
def assembly(e):
    # Get element stiffness and force
    ke, fe = heat2delem(e)

    LM = setup_ID_LM (var.neq, var.nel, var.nd)
    for loop1 in range(var.nen):
        i = LM[loop1][e]
        var.f[i] = var.f[i] + fe[loop1] # Assemble forces
        for loop2 in range(var.nen):
            j = LM[loop2][e]
            var.K[i][j] = var.K[i][j] + ke[loop1][loop2] # Assemble stiffness  
    return

#=============================================================
#               source and flux (src_and_flux.m)
#=============================================================
def src_flux(neq, nbe, ngp):
    for i in range(0,neq):
        var.f[var.ID[i]] = var.f[var.ID[i]] + var.P[var.ID[i]]

    x, y = physCoord(var.lpx, var.lpy)
    for i in range(0,nbe):
        fq = np.zeros((2,1))
        n_bce = np.zeros((2,1))

        node1 = int(var.n_bc[0][i])
        node2 = int(var.n_bc[0][i+1])
        n_bce[0] = var.n_bc[1][i]
        n_bce[1] = var.n_bc[1][i+1]

        x1 = x[node1]
        y1 = y[node1]
        x2 = x[node2]
        y2 = y[node2]

        length = math.sqrt((x2-x1)**2 + (y2-y1)**2)
        J = length/2

        w,gp = gauss(ngp)

        for i in range(0,ngp):
            N = np.zeros((2,1))
            psi = gp[i]
            N[0] = 0.5*(1-psi)
            N[1] = 0.5*(1+psi)

            flux = N[0]*n_bce[0] + N[1]*n_bce[1]
            fq[0] = fq[0] + w[i]*N[0]*flux*J
            fq[1] = fq[1] + w[i]*N[1]*flux*J

        fq = -fq

        var.f[var.ID[node1]] = var.f[var.ID[node1]] + fq[0]
        var.f[var.ID[node2]] = var.f[var.ID[node2]] + fq[1]
    return 



#=============================================================
#               Solve the system (solvedr.m)
#=============================================================




#=============================================================
#               get flux (get_flux.m)
#=============================================================




#=============================================================
#               Postprocess (postprocessor.m)
#=============================================================