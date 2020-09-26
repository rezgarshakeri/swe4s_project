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
nelx = 3                      # number of elements in x direction
nely = 2                      # number of elements in y direction
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

ID   = np.zeros((neq,1))
d    = np.zeros((neq,1))
# gauss Integration
ngp  = 2                      # number of gauss points

# boundary conditions on left and bottom
nd = lpx+lpy -1               # number of nodes on essential boundary

# essential B.C. (prescribed temperature)
for i in range(0,lpx):
    flags [i] = 2
    e_bc  [i] = 10.0        # bottom edge

for i in range(lpx,nnp-nelx,lpx):
    flags[i] = 2
    e_bc[i]  = -10.0         # left edges

# natural B.C  - defined on edges positioned on natural boundary
n_bc = np.zeros((4, nelx))
nbc = nnp-nelx
for i in range(0,nelx):
    n_bc[0][i] = nbc + i
    n_bc[1][i] = nbc + 1 + i
    n_bc[2][i] = 0
    n_bc[3][i] = 0

nbe = nelx

