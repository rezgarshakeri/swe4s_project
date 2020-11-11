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
    # we test the connectivity matrix for nelx=2, nely=1, as above
    # connectivity returns an array which relate the local
    # numbering to global numbering. Thenumber of rows of connectivity
    # is 4 since we use 4-nodes element and number of columns of connectivity
    # is equal to number of element we have.
    # for example compare the node numbering of element (e) with
    # element (1): 0-->1, 1-->2, 2-->5, 3-->4
    # so second column of connectivity will be 1,2,5,4
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
            self.assertEqual(flags[i][0], 2)


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


def main():
    unittest.main()


if __name__ == '__main__':
    main()
