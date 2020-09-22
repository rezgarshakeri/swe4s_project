import numpy as np

#=============================================================
#               Materials properties and mesh
#=============================================================

# material property
k  = 5                        # thermal conductivity
D  = k*np.identity(2)         # conductivity matrix
    
# mesh specifications    
nsd  = 2                      # number of space dimensions
nen  = 4                      # number of element nodes 
ndof = 1                      # degrees-of-freedom per node
nelx = 2                      # number of elements in x direction
nely = 3                      # number of elements in y direction
nel  = nelx*nely              # number of elements
lpx = nelx+1                  # number of nodes in the x direction 
lpy = nely+1                  # number of nodes in the y direction 
nnp  = (lpx)*(lpy)            # number of nodal nodes    
neq  = nnp*ndof               # number of equations

f = np.zeros((neq,1))         # initialize nodal source vector
d = np.zeros((neq,1))         # initialize nodal temperature vector
K = np.zeros((neq,neq))       # initialize stiffness matrix

flags= np.zeros((neq,1))      # array to set B.C flags 
e_bc = np.zeros((neq,1))      # essential B.C array
P    = np.zeros((neq,1))      # initialize point source defined at a node
s    = 6*np.ones((nen,nel))   # heat source

# gauss Integration
ngp  = 2                      # number of gauss points

# boundary conditions on left and bottom
nd = lpx+lpy -1               # number of nodes on essential boundary

# essential B.C. (prescribed temperature)
for i in range(0,lpx):
    flags [i,0] = 2
    e_bc  [i,0] = 0.0        # bottom edge

for i in range(lpx,nnp-nelx,lpx):
    flags[i,0] = 2
    e_bc[i,0]  = 0.0         # left edges

# natural B.C  - defined on edges positioned on natural boundary
n_bc = np.zeros((4, nelx))
nbc = nnp-nelx
for i in range(0,nelx):
    n_bc[0,i] = nbc + i
    n_bc[1,i] = nbc + 1 + i
    n_bc[2,i] = 20
    n_bc[3,i] = 20

nbe = nelx

#=============================================================
#               Mesh Generation and IEN (mesh2d.m)
#=============================================================

# This function returns the physical coordinates of the nodes
def physCoord(lpx, lpy):
    
    nnp = lpx * lpy
    x0 = np.linspace(0,1,lpx)
    y0 = 0.5*x0**2                    # the bottom geometry line

    y = np.zeros((nnp,1))
    for i in range(0,lpx):
        y1 = np.linspace(y0[i],1,lpy)
        for j in range(0,lpy):        
            y[i + j*lpx] = y1[j]     # collection of y coordinate

    x = np.zeros((nnp,1))
    for i in range(0,lpy):        
        for j in range(0,lpx):
            x[j + i*lpx] = x0[j]   # collection of x coordinate
    return x, y


# generate the IEN connectivity array
def connectivity(nel, lpx):

    IEN = np.zeros((4,nel))
    rowcount = 0
    for elementcount in range(1,nel+1):
        IEN[0,elementcount-1] = elementcount + rowcount
        IEN[1,elementcount-1] = elementcount + 1 + rowcount
        IEN[2,elementcount-1] = elementcount + (lpx + 1) + rowcount
        IEN[3,elementcount-1] = elementcount + (lpx) + rowcount
        if np.mod(elementcount,lpx-1) == 0:
            rowcount = rowcount + 1
            
    return IEN

#=============================================================
#               Plot Mesh (plotmesh.m)
#=============================================================

 


#=============================================================
#               Setup ID and LM (setup_ID_LM.m)
#=============================================================




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

    detJ = J[0,0]*J[1,1] - J[0,1]*J[1,0]
    invJ = np.zeros((2,2))
    invJ[0,0] =  J[1,1]
    invJ[0,1] = -J[0,1]
    invJ[1,0] = -J[1,0]
    invJ[1,1] =  J[0,0]
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
    return


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