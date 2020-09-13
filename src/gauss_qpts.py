
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