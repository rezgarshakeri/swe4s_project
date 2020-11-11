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
    def test_phys_coord_4element(self):
        # at x=0.5 ==> y=0.5*(0.5)^2=0.125
        node1_y = 0.125
        node4_x = 0.5
        x, y = FE.phys_coord(2, 2)
        # the output is x (9x1 array) and y (9x1 array)
        self.assertEqual(y[1][0], node1_y)
        self.assertEqual(x[4][0], node4_x)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
