import numpy as np
# material property
k  = 5;        # thermal conductivity
D  = k*np.identity(2); # conductivity matrix

# mesh specifications
nsd  = 2;         # number of space dimensions
nen  = 4;         # number of element nodes 
ndof = 1;         # degrees-of-freedom per node
nelx = 2;         # number of elements in x direction
nely = 2;         # number of elements in y direction
nel  = nelx*nely;         # number of elements
lpx = nelx+1;             # number of nodes in the x direction 
lpy = nely+1;             # number of nodes in the y direction 
nnp  = (lpx)*(lpy);       # number of nodal nodes

neq  = nnp*ndof;  # number of equations


f = np.zeros((neq,1));      # initialize nodal source vector
d = np.zeros((neq,1));      # initialize nodal temperature vector
K = np.zeros((neq,neq));        # initialize stiffness matrix

flags = np.zeros((neq,1));  # array to set B.C flags 
e_bc  = np.zeros((neq,1));  # essential B.C array
#n_bc  = np.zeros((neq,1));  # natural B.C array
P    = np.zeros((neq,1));   # initialize point source defined at a node
s    = 6*np.ones((nen,nel));  # heat source

# gauss Integration
ngp    = 2;                          # number of gauss points

# boundary conditions and point forces
nd = lpx+nely;     # number of nodes on essential boundary

# essential B.C.
flags(1:lpx)    = 2;                e_bc(1:lpx)     = 0.0;  # bottom edge
flags(lpx+1:lpx:nnp-nelx) = 2;      e_bc(lpx+1:lpx:nnp-nelx)  = 0.0;  #left edge


# plots
#compute_flux = 'yes';
#plot_mesh    = 'yes';
#plot_nod     = 'yes';
#plot_temp    = 'yes';
#plot_flux    = 'yes';

# natural B.C  - defined on edges positioned on natural boundary
n_bc = np.zeros((4, nelx))
n_bc(1,1:nelx) = nnp-nelx:nnp-1
n_bc(2,1:nelx) = nnp-nelx+1:nnp
n_bc(3,1:nelx) = 20
n_bc(4,1:nelx) = 20

nbe = nelx


# mesh generation
# mesh2d;


 
 