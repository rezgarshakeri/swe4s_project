import numpy as np

# N defines bilinear basis function
# N = [N1, N2, N3, N4] 
def basis(xi, eta):
    N = 0.25 * np.array([(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)])
    return N