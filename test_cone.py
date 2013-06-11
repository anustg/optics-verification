# Test that the infinite and finite cone geometry does what it should.

import unittest
import numpy as N
N.set_printoptions(linewidth=150)

from tracer.ray_bundle import RayBundle
from tracer.cone import InfiniteCone, Cone
from tracer.spatial_geometry import rotx

class TestInfiniteCone(unittest.TestCase):
    def setUp(self):
        pos = N.c_[[-1.,0,0], [2,0,0], [1,0,2]]
        dire = N.c_[[0.,0,1], [-1,0,1], [-1,2,0]]
        self.prm = N.r_[1,N.sqrt(2),N.sqrt(5)]
        self.pts = N.c_[[-1.,0,1],[1,0,1],[0,2,2]]
        self.nrm = N.c_[[-1.,0,-1],[1,0,-1],[0,-1,1]]

        dire /= N.sqrt(N.sum(dire**2, axis=0)) # normalise dire
        self.nrm /= N.sqrt(N.sum(self.nrm**2, axis=0)) # normalise nrm

        self.bund = RayBundle(vertices=pos, directions=dire)
        self.gm = InfiniteCone(c = 1.)

    def test_as_placed1(self):

        prm = self.gm.find_intersections(N.eye(4), self.bund)
        #print "prm",prm
        N.testing.assert_array_almost_equal(prm, self.prm)

        self.gm.select_rays(N.arange(0,3))

        nrm = self.gm.get_normals()
        #print "self.nrm\n",self.nrm
        #print "nrm\n",nrm
        N.testing.assert_array_almost_equal(nrm, self.nrm)

        pts = self.gm.get_intersection_points_global()

        print "pts",pts
        N.testing.assert_array_almost_equal(pts, self.pts)


class TestInfiniteConeShifted(unittest.TestCase):
    """
    Very similar to TestInfiniteCone, but add 'a=-1' to shift the cone
    apex to z=-1 on the z axis, check that everything still works correctly.
    """
    def setUp(self):
        pos = N.c_[[-1.,0,-1], [2,0,-1], [1,0,1]]
        dir = N.c_[[0.,0,1], [-1,0,1], [-1,2,0]]
        self.prm = N.r_[1,N.sqrt(2),N.sqrt(5)]
        self.pts = N.c_[[-1.,0,0],[1,0,0],[0,2,1]]
        self.nrm = N.c_[[-1.,0,-1],[1,0,-1],[0,-1,1]]

        dir /= N.sqrt(N.sum(dir**2, axis=0)) # normalise dir
        self.nrm /= N.sqrt(N.sum(self.nrm**2, axis=0)) # normalise nrm

        self.bund = RayBundle(vertices=pos, directions=dir)
        self.gm = InfiniteCone(c = 1., a = -1.)

    def test_as_placed2(self):

        prm = self.gm.find_intersections(N.eye(4), self.bund)
        #print "prm",prm
        N.testing.assert_array_almost_equal(prm, self.prm)

        self.gm.select_rays(N.arange(0,3))

        nrm = self.gm.get_normals()
        #print "self.nrm\n",self.nrm
        #print "nrm\n",nrm
        N.testing.assert_array_almost_equal(nrm, self.nrm)

        pts = self.gm.get_intersection_points_global()

        print "pts",pts
        N.testing.assert_array_almost_equal(pts, self.pts)


class TestCone(unittest.TestCase):
    def setUp(self):
        pos = N.c_[[-1.,0,0], [2,0,0], [1,0,2], [4,0,0],[2,0,1]]
        dir = N.c_[[0.,0,1], [-1,0,1], [-1,2,0], [0,0,1],[-1,0,0]]
        self.prm = N.r_[1,N.sqrt(2),N.sqrt(5),N.inf,1]
        self.pts = N.c_[[-1.,0,1],[1,0,1],[0,2,2],[1,0,1]]
        self.nrm = N.c_[[-1.,0,-1],[1,0,-1],[0,-1,1],[1,0,-1]]

        dir /= N.sqrt(N.sum(dir**2, axis=0)) # normalise dir
        self.nrm /= N.sqrt(N.sum(self.nrm**2, axis=0)) # normalise nrm

        self.bund = RayBundle(vertices=pos, directions=dir)
        self.gm = Cone(r = 3., h = 3.)

    def test_as_placed3(self):

        #print "ASKING FOR INTERSECTIONS\n"
        prm = self.gm.find_intersections(N.eye(4), self.bund)
        #print "GOT INTERSECTIONS\n"
        #print "prm=\n",prm
        N.testing.assert_array_almost_equal(prm, self.prm)

        self.gm.select_rays(self.prm != N.inf)

        nrm = self.gm.get_normals()
        #print "self.nrm\n",self.nrm
        #print "nrm\n",nrm
        N.testing.assert_array_almost_equal(nrm, self.nrm)

        pts = self.gm.get_intersection_points_global()

        print "pts",pts
        N.testing.assert_array_almost_equal(pts, self.pts)



if __name__ == '__main__':
    unittest.main()


# vim: et:ts=4
