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
for i1 in range(0,lpx):
    flags [i1,0] = 2
    e_bc  [i1,0] = 0.0        # bottom edge

for i2 in range(lpx,nnp-nelx,lpx):
    flags[i2,0] = 2
    e_bc[i2,0]  = 0.0         # left edges

# natural B.C  - defined on edges positioned on natural boundary
n_bc = np.zeros((4, nelx))
nbc = nnp-nelx
for i3 in range(0,nelx):
    n_bc[0,i3] = nbc + i3
    n_bc[1,i3] = nbc + 1 + i3
    n_bc[2,i3] = 20
    n_bc[3,i3] = 20

nbe = nelx

#=============================================================
#               Mesh Generation and IEN (mesh2d.m)
#=============================================================

x0 = np.linspace(0,1,lpx)
y0 = 0.5*x0**2                    # the bottom geometry line

y = np.zeros((nnp,1))
for i4 in range(0,lpx):
    y1 = np.linspace(y0[i4],1,lpy)
    for j in range(0,lpy):        
        y[i4 + j*lpx] = y1[j]     # collection of y coordinate


x = np.zeros((nnp,1))
for i5 in range(0,lpy):        
    for j1 in range(0,lpx):
        x[j1 + i5*lpx] = x0[j1]   # collection of x coordinate


# generate the IEN connectivity array
IEN = np.zeros((4,nel))
rowcount = 0
for elementcount in range(1,nel+1):
    IEN[0,elementcount-1] = elementcount + rowcount
    IEN[1,elementcount-1] = elementcount + 1 + rowcount
    IEN[2,elementcount-1] = elementcount + (lpx + 1) + rowcount
    IEN[3,elementcount-1] = elementcount + (lpx) + rowcount
    if np.mod(elementcount,lpx-1) == 0:
        rowcount = rowcount + 1

#=============================================================
#               Plot Mesh (plotmesh.m)
#=============================================================

 


#=============================================================
#               Setup ID and LM (setup_ID_LM.m)
#=============================================================




#=============================================================
#               basis and derivative of basis and gauss
#=============================================================

# basis here

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

# gauss here 


#=============================================================
#               stiffness, force element (heat2Delem.m)
#=============================================================





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