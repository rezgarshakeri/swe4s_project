import numpy as np
# derivativ of basis function, 
# [[N1_xi,N2_xi,N3_xi,N4_xi],
#  [N1_eta,N2_eta,N3_eta,N4_eta]]
# coord is the coordinate of the element's nodes in physical coordinate
# coord = [[x1,y1],[x2,y2],[x3,y3],[x4,y4]]; (4x2)
# topic3,pg.11: x(xi,eta)=[xi,eta]^T, xe=[x,y]^T, x = N(xi,eta)*xe
# topic3,pg.13: J = dx/d(xi,eta) = dN/d(xi,eta) *xe
#finally B = dN/dx = dN/d(xi,eta) * d(xi,eta)/dx = inv(J)*dN/d(xi,eta)
def d_basis(xi,eta,coord):

    #Calculate the Grad(N) matrix
    #[[],
    #  []]
    dN = np.array(0.25*[[eta-1, 1-eta, 1+eta,-eta-1],
                        [xi-1 , -xi-1, 1+xi , 1-xi]])

    J     = dN*coord      # compute Jacobian matrix 
    detJ  = np.linalg.det(J)     # Jacobian
      
    B     = np.linalg.solve(J, dN)       # compute the B matrix

    return B, detJ