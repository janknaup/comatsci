'''
Created on Mar 13, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
from __future__ import print_function
import unittest,numpy,os
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
        g.writexyz("test.xyz")
        g2=Geometry.Geometry()
        g2.readfile("test.xyz")
        os.unlink("test.xyz")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
    
    
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
        g.writegen("test.gen")
        g2=Geometry.Geometry()
        g2.readfile("test.gen")
        os.unlink("test.gen")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        
    
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
        g.writeCDH("test.cdh")
        g2=Geometry.Geometry()
        g2.readCDHFile("test.cdh")
        os.unlink("test.cdh")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        
    
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
        g.writefmg("test.fmg")
        g2=Geometry.Geometry()
        g2.readfile("test.fmg")
        os.unlink("test.fmg")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        
        
    def testAimsOut(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        g.writeAIMS("test.aims")
        testfile=open("test.aims")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.aims")
        self.assertEqual(g.getAimsString(), "".join(testlines).strip("\n"), "written file corrupted")
        AtomTypes=["C","C","H","H","H","H"]
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        georef*=Geometry.Angstrom
        for ii in range(g.Atomcount):
            tokens=testlines[ii+1].split()
            self.assertEqual(AtomTypes[ii],tokens[4],"Element type differs")
            self.assertAlmostEqual(georef[ii][0],float(tokens[1]),6,"X coordinate differs")
            self.assertAlmostEqual(georef[ii][1],float(tokens[2]),6,"Y coordinate differs")
            self.assertAlmostEqual(georef[ii][2],float(tokens[3]),6,"Z coordinate differs")


    def testTmOut(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        g.writeTurboMole("test.tm")
        testfile=open("test.tm")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.tm")
        self.assertEqual(g.tmString().strip(), "".join(testlines).strip(), "written file corrupted")
        AtomTypes=["C","C","H","H","H","H"]
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        for ii in range(g.Atomcount):
            tokens=testlines[ii+1].split()
            self.assertEqual(AtomTypes[ii],tokens[3],"Element type differs")
            self.assertAlmostEqual(georef[ii][0],float(tokens[0]),6,"X coordinate differs")
            self.assertAlmostEqual(georef[ii][1],float(tokens[1]),6,"Y coordinate differs")
            self.assertAlmostEqual(georef[ii][2],float(tokens[2]),6,"Z coordinate differs")
    
    
    def testFdfOut(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        g.writefdf("test.fdf")
        testfile=open("test.fdf")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.fdf")
        reffile=open("../../../test_regr/unittest_data/geotest1.fdf")
        reflines=list(reffile)
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written fdf file differs from reference")


    def testPDBOut(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        g.writepdb("test.pdb")
        testfile=open("test.pdb")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.pdb")
        reffile=open("../../../test_regr/unittest_data/geotest1.pdb")
        reflines=list(reffile)
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written pdb file differs from reference")
            
            
    def testXYZQOut(self):
        g=Geometry.Geometry()
        g.readfile("../../../test_regr/unittest_data/geotest1.fmg")
        g.AtomCharges=[0.,1.,3.,7.8220452,5.5174032,-3.14157324839284729457398775532]
        g.writexyzq("test.xyzq")
        testfile=open("test.xyzq")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.xyzq")
        self.assertEqual(g.pointchargesstring().strip(), "".join(testlines).strip(), "written file corrupted")
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        for ii in range(g.Atomcount):
            tokens=testlines[ii].split()
            self.assertAlmostEqual(georef[ii][0],float(tokens[0]),6,"X coordinate differs")
            self.assertAlmostEqual(georef[ii][1],float(tokens[1]),6,"Y coordinate differs")
            self.assertAlmostEqual(georef[ii][2],float(tokens[2]),6,"Z coordinate differs")
            self.assertAlmostEqual(g.AtomCharges[ii],float(tokens[3]),6,"charge differs")


def suite():
    suite = unittest.makeSuite(geometryTest,'test')
    return suite


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testGeometry']
    unittest.main()