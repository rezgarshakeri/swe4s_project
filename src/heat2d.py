""" This code solves 2d steady state heat equation
using finite element numerical method.
The output is the nodal Temperature.
"""

import FE_subroutines as FE
import plot
import argparse


def main():
    parser = argparse.ArgumentParser(description='Input for FE_subroutines')

    parser.add_argument('--num_elm_x',
                        dest='num_elm_x',
                        type=int,
                        default=2,
                        help='Number of element in x direction')

    parser.add_argument('--num_elm_y',
                        dest='num_elm_y',
                        type=int,
                        default=3,
                        help='Number of element in y direction')

    parser.add_argument('--T0_bottom',
                        dest='T0_bottom',
                        type=float,
                        default=10,
                        help='Prescribed Temperature on the bottom edge')

    parser.add_argument('--T0_left',
                        dest='T0_left',
                        type=float,
                        default=-10,
                        help='Prescribed Temperature on the left edge')

    parser.add_argument('--heat_source',
                        dest='heat_source',
                        type=float,
                        default=6,
                        help='Heat source on the domain')

    parser.add_argument('--flux_top',
                        dest='flux_top',
                        type=float,
                        default=0,
                        help='Prescribed flux on the top edge')
    args = parser.parse_args()

    nelx = args.num_elm_x
    nely = args.num_elm_y
    T0_bottom = args.T0_bottom
    T0_left = args.T0_left
    s0 = args.heat_source
    flux_top = args.flux_top

    # Gauss point, 2 is better than 1 since it has less errors
    ngp = 2

    # Assemble stiffness and forcing vector
    K, f = FE.assembly(nelx, nely, ngp, s0, T0_bottom, T0_left)

    # Update forcing vector f by adding flux vector
    F = FE.src_flux(nelx, nely, T0_bottom, T0_left, flux_top, ngp, f)

    # Get the prescribed Temperature on the boundary nodes
    d0, ID, LM = FE.setup_ID_LM(nelx, nely, T0_bottom, T0_left)

    # Find the nodal Temperature on the domain
    d = FE.solvedr(nelx, nely, K, F, d0)

    # Print the nodal Temperature
    print(d)

    # Plot the 2d temperature distribution
    plot.plot_temp(nelx, nely, T0_bottom, T0_left, K, F)

    return


if __name__ == '__main__':
    main()
