import unittest
import FE_subroutines as FE


class Testphys_coord(unittest.TestCase):
    # Consider only 1 element (nelx = 1, nely = 1).
    #  So we have 4 nodes, start from the origin of x-y coordinate,
    # and y=0.5x^2 for bottom edge (see documention of phy_coord function)
    # coordinate: node0=(0,0), node1=(1,0.5), node2=(0,1), node3=(1,1)
    def test_phys_coord_1element(self):
        node0_x = 0.
        node1_y = 0.5
        x, y = FE.phys_coord(1, 1)
        # the output is x (4x1 array) and y (4x1 array)
        self.assertEqual(x[0][0], node0_x)
        self.assertEqual(y[1][0], node1_y)

    # Now we test for 4 elements (nelx=2, nely=2), exactly like
    # the figure in documentation of phys_coord
    def test_phys_coord_4elements(self):
        # at x=0.5 ==> y=0.5*(0.5)^2=0.125
        node1_y = 0.125
        node4_x = 0.5
        x, y = FE.phys_coord(2, 2)
        # x (9x1 array) and y (9x1 array) since we have 9 nodes
        self.assertEqual(y[1][0], node1_y)
        self.assertEqual(x[4][0], node4_x)


class Testconnectivity(unittest.TestCase):
    """
    (for nelx = 2, nely = 1)
    Global numbering of nodes
    3---------4----------5
    |         |   (1)    |
    |   (0)   |      ----2
    |      ---1-----/
    0-----/
    Local numbering for element (e)
    3-----2
    | (e) |
    0-----1
    """
    # We test the connectivity matrix for nelx=2, nely=1, as above.
    # Connectivity returns an array which relates the local
    # numbering to the global numbering. The number of rows of the connectivity
    # is 4 since we use 4-node elements and the number of the columns of
    # the connectivity is equal to number of elements we have.
    # For example, compare the node numbering of element (e) with
    # element (1): 0-->1, 1-->2, 2-->5, 3-->4
    # So the second column of the connectivity will be 1,2,5,4
    def test_connectivity_2elements(self):
        IEN_test = [[0, 1], [1, 2], [4, 5], [3, 4]]
        IEN = FE.connectivity(2, 1)

        for i in range(0, 4):
            for j in range(0, 2):
                self.assertEqual(IEN_test[i][j], IEN[i][j])


class TestDirichlet_BCs(unittest.TestCase):
    # we test with same figure as we have for connectivity
    # so, we set the Dirichlet boundary condition of bottom
    # edge to 10 and left edge to -10.
    # so we should have 3 nodes (0,1,2) with 10
    # and 1 node (3) with -10 in e_bc array
    # flags must be [2,2,2,2,0,0], since we have 4 nodes with BCs
    def test_Dirichlet_BCs_2elements(self):
        test_bottom = 10
        test_left = -10
        e_bc, flags = FE.Dirichlet_BCs(2, 1, 10, -10)
        # test bottom edge
        for i in range(0, 3):
            self.assertEqual(e_bc[i][0], test_bottom)
        # test left edge (for the above mesh we have only 1 nodes)
        self.assertEqual(e_bc[3][0], test_left)

        # test flags
        for i in range(0, 4):
            self.assertEqual(flags[i][0], 1)


class Testsetup_ID_LM(unittest.TestCase):
    # see the documentation and code for what setup_ID_LM should returns
    def test_setup_ID_LM_2elements(self):
        d0_test = [[10.], [10.], [10.], [-10.], [0.], [0.]]
        ID_test = [[0], [1], [2], [3], [4], [5]]
        LM_test = [[0, 1], [1, 2], [4, 5], [3, 4]]
        d0, ID, LM = FE.setup_ID_LM(2, 1, 10, -10)
        # test initial temperature in the domain
        for i in range(0, 6):
            self.assertEqual(d0[i][0], d0_test[i][0])
        # test ID array
        for i in range(0, 6):
            self.assertEqual(ID[i][0], ID_test[i][0])
        # test LM array
        for i in range(0, 4):
            for j in range(0, 2):
                self.assertEqual(LM[i][j], LM_test[i][j])


class Testbasis(unittest.TestCase):
    """
    3-----2
    | (e) |
    0-----1
    """
    # basis return and 4x1 array as N = [N0 N1 N2 N3], for nodes 0 to 3
    # Origin of (xi, eta) coordinate is
    # at the center of element e. the coordinates of the 4 nodes are
    # node0 = (-1, 1), node1 = (1, -1), node2 = (1, 1), node3 = (-1, 1)
    # for N = basis(-1, 1) N will be N = [[1. 0. 0. 0.]]
    # for N = basis(1, 1) N will be N = [[0. 0. 1. 0.]]
    # test shape function of node0
    def test_basis0(self):
        N_test0 = [[1., 0., 0., 0.]]
        N = FE.basis(-1, -1)

        for i in range(0, 4):
            self.assertEqual(N[0][i], N_test0[0][i])

    # test shape function of node1
    def test_basis1(self):
        N_test1 = [[0., 1., 0., 0.]]
        N = FE.basis(1, -1)

        for i in range(0, 4):
            self.assertEqual(N[0][i], N_test1[0][i])

    # test shape function of node2
    def test_basis2(self):
        N_test2 = [[0., 0., 1., 0.]]
        N = FE.basis(1, 1)

        for i in range(0, 4):
            self.assertEqual(N[0][i], N_test2[0][i])

    # test shape function of node3
    def test_basis3(self):
        N_test3 = [[0., 0., 0., 1.]]
        N = FE.basis(-1, 1)

        for i in range(0, 4):
            self.assertEqual(N[0][i], N_test3[0][i])


class Testheat2d(unittest.TestCase):
    # test the final solution
    def test_heat2d_1(self):

        nelx = 2
        nely = 3
        T0_bottom = 10
        T0_left = -10
        flux_top = 0
        ngp = 2
        s0 = 6
        # assemble stiffness and forcing vector
        K, f = FE.assembly(nelx, nely, ngp, s0, T0_bottom, T0_left)
        # update forcing vector f by adding flux vector
        F = FE.src_flux(nelx, nely, T0_bottom, T0_left, flux_top, ngp, f)
        # get the prescribed temperature on the boundary nodes
        d0, ID, LM = FE.setup_ID_LM(nelx, nely, T0_bottom, T0_left)
        # find the nodal temperature on the domain
        d = FE.solvedr(nelx, nely, K, F, d0)
        # since we have nelx = 2, the first 3 entries of the d must be
        # equal to prescribed temperature of the bottom edge
        for i in range(0, nelx + 1):
            self.assertEqual(d[i][0], T0_bottom)

    def test_heat2d_2(self):

        nelx = 2
        nely = 2
        T0_bottom = 5
        T0_left = -2
        flux_top = 0
        ngp = 2
        s0 = 0
        # assemble stiffness and forcing vector
        K, f = FE.assembly(nelx, nely, ngp, s0, T0_bottom, T0_left)
        # update forcing vector f by adding flux vector
        F = FE.src_flux(nelx, nely, T0_bottom, T0_left, flux_top, ngp, f)
        # get the prescribed temperature on the boundary nodes
        d0, ID, LM = FE.setup_ID_LM(nelx, nely, T0_bottom, T0_left)
        # find the nodal temperature on the domain
        d = FE.solvedr(nelx, nely, K, F, d0)
        # since we have nelx = 2, the first 3 entries of the d must be
        # equal to prescribed temperature of the bottom edge
        for i in range(0, nelx + 1):
            self.assertEqual(d[i][0], T0_bottom)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
