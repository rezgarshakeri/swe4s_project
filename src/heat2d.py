import math
import numpy as np

import variables as var
import FE_subroutines as FE


def main():

    # Forms the assembly matrix by looping over all elements
    for e in range(var.nel):
        FE.assembly(e)

    # Adds the forcing term causing by heat flux
    FE.src_flux(var.neq, var.nbe, var.ngp)

    # Solves the linear system
    d = FE.solvedr(var.neq, var.nd)

    # Prints out the Temperatures at each node
    # TODO: We need to create plots instead of printing
    print(d)

    return


if __name__ == '__main__':
    main()
