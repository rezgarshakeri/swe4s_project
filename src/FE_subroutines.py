import numpy as np
import variables as vr

#=============================================================
#               Mesh Generation and IEN (mesh2d.m)
#=============================================================

# This function returns the physical coordinates of the nodes
def physCoord(lpx, lpy):

    x0 = np.linspace(0,1,lpx)
    y0 = 0.5*x0**2                    # the bottom geometry line

    y = np.zeros((vr.nnp,1))
    for i in range(0,lpx):
        y1 = np.linspace(y0[i],1,lpy)
        for j in range(0,lpy):        
            y[i + j*lpx] = y1[j]     # collection of y coordinate

    x = np.zeros((vr.nnp,1))
    for i in range(0,lpy):        
        for j in range(0,lpx):
            x[j + i*lpx] = x0[j]   # collection of x coordinate
    return x, y


# generate the IEN connectivity array
def connectivity(nel, lpx):

    IEN = np.zeros((4,nel))
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
def setup_ID_LM (neq,nel,nd):

    IEN = connectivity(nel,vr.lpx)
    LM = np.zeros((4,nel))
    count  = 0
    count1 = 0
    for i in range(0,neq):
        if vr.flags[i] == 2:
            vr.ID[i]  = count
            vr.d[count] = vr.e_bc[i]
            count  = count +1
        else:
            vr.ID[i]  = count1 + nd
            count1 = count1 + 1
  
    for i in range(0,nel):
        for j in range(0,vr.nen):
            ien = int(IEN[j][i])
            LM[j][i] = vr.ID[ien]

    return LM

#=============================================================
#               basis and derivative of basis and gauss
#=============================================================

# N defines bilinear basis function
# N = [N1, N2, N3, N4] 
def basis(xi, eta):
    N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
    return N


def d_basis(xi,eta,coord):

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
        return gp,w
    elif ngp == 2:
        gp = [-0.57735027,  0.57735027 ]  
        w  = [1,            1          ]
        return gp,w
    else:
        print("Error: This code supports only 1 and 2 quadrature points.")
    


#=============================================================
#               stiffness, force element (heat2Delem.m)
#=============================================================
# TODO
# e: number of element
def heat2delem(e):

    ke = np.zeros(nen, nen) # Initialize element conductance matrix
    fe = np.zeros(nen, 1)   # Initialize element nodal source vector

    # Get coordinates of element nodes 
        # TODO: define IEN, x, and y in inputData
    je = IEN(:,e)                       
    C  = np.transpose( [ [x(je)] , [y(je)] ] )

    # Get gauss points and weights
    w, gp = gauss(nd.ngp)




#=============================================================
#               Assembly (assembly.m)
#=============================================================




#=============================================================
#               source and flux (src_and_flux.m)
#=============================================================




#=============================================================
#               Solve the system (solvedr.m)
#=============================================================




#=============================================================
#               get flux (get_flux.m)
#=============================================================




#=============================================================
#               Postprocess (postprocessor.m)
#=============================================================