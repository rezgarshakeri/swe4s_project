# e: number of element

import inputData as nd
import numpy as np
import gauss_qpts as gp

def heat2delem(e):

    ke = np.zeros(nen, nen) # Initialize element conductance matrix
    fe = np.zeros(nen, 1)   # Initialize element nodal source vector

    # Get coordinates of element nodes 
        # TODO: define IEN, x, and y in inputData
    je = nd.IEN(:,e)                       
    C  = np.transpose( [ [nd.x(je)] , [nd.y(je)] ] )

    # Get gauss points and weights
    w, gp = gp.gauss(nd.ngp)

global ndof nnp nel nen nsd neq ngp nee neq 
global nd e_bc s P D
global LM ID IEN flags n_bc 
global x y nbe
global compute_flux plot_mesh plot_temp plot_flux plot_nod
global nelx nely lpx lpy