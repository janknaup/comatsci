'''
Created on Mar 13, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
import unittest,numpy
from comatsci import Geometry

class geometryTest(unittest.TestCase):


    def testGeometryBasics(self):
        g=Geometry.Geometry()
        self.assertTrue(isinstance(g, Geometry.Geometry), "Geometry Constructor does not return Geometry instance")
        g.addatom(1,(0.,0.,0.))
        g.addatom(8,(2,0.,0.))
        self.assertEqual(g.AtomTypes, [1,8], "Atom Types differ")
        georef=numpy.array([[0.,0.,0.],[2.,0.,0.]],dtype=float)
        self.assertTrue(numpy.array_equal(g.Geometry, georef), "Coordinates array differs")
        self.assertTrue(numpy.array_equal(g.AtomCharges,[0.,0.]), "Charges array differs")
        self.assertEqual(["H","O"], g.AtomSubTypes, "geometry subtypes differ")
        g.delatom(0)
        georef=numpy.array([[2.,0.,0.],],dtype=float)
        self.assertEqual(g.AtomTypes, [8,], "After delatom: types differ")
        self.assertTrue(numpy.array_equal(g.Geometry, georef), "After delatom: coordinates differ")
        g.delatom(0)
        ## TODO: ADD GeometryError when trying to delete nonexistent atom
        ## self.assertRaises(Geometry.GeometryError, g.delatom(0))
        waterlayer=g.addlayer("water")
        self.assertEqual("water", g.LayerDict[waterlayer].Name, "added layer name differs")
        g.addatom(1,(0.,0.,0.),waterlayer)
        g.addatom(8,(2.,0.,0.),waterlayer)
        lsg=g.layersubgeometry(waterlayer)
        self.assertEqual(g.AtomTypes, lsg.AtomTypes, "layersubgeometry types differ")
        self.assertTrue(numpy.array_equal(g.Geometry, lsg.Geometry), "layersubgeometry coordinates differ")
        self.assertEqual(g.AtomCharges, lsg.AtomCharges, "layersubgeometry charges differ")
        
    
    def testXYZio(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.xyz")
        self.assertEqual(g.Mode, "C", "geometry mode differs")
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g.AtomTypes, "geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g.AtomSubTypes, "geometry subtypes differ")

    
    
    def testGENio(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.gen")
        self.assertEqual(g.Mode, "C", "geometry mode differs")
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g.AtomTypes, "geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g.AtomSubTypes, "geometry subtypes differ")
        
    
    def testCDHio(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.gen")
        self.assertEqual(g.Mode, "C", "geometry mode differs")
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g.AtomTypes, "geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g.AtomSubTypes, "geometry subtypes differ")
        
    
    def testFMGio(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        self.assertEqual(g.Mode, "C", "geometry mode differs")
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g.AtomTypes, "geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g.AtomSubTypes, "geometry subtypes differ")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testGeometry']
    unittest.main(verbosity=2)