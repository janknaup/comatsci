'''
Created on Mar 13, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
from __future__ import print_function
import unittest,numpy,os
import comatsci

class geometryTest(unittest.TestCase):


    def testGeometryBasics(self):
        """basic Geometry functionality"""
        g=comatsci.geometry.Geometry()
        self.assertTrue(isinstance(g, comatsci.geometry.Geometry), "Geometry Constructor does not return Geometry instance")
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
        ## self.assertRaises(comatsci.geometry.GeometryError, g.delatom(0))
        waterlayer=g.addlayer("water")
        self.assertEqual("water", g.LayerDict[waterlayer].Name, "added layer name differs")
        g.addatom(1,(0.,0.,0.),waterlayer)
        g.addatom(8,(2.,0.,0.),waterlayer)
        lsg=g.layersubgeometry(waterlayer)
        self.assertEqual(g.AtomTypes, lsg.AtomTypes, "layersubgeometry types differ")
        self.assertTrue(numpy.array_equal(g.Geometry, lsg.Geometry), "layersubgeometry coordinates differ")
        self.assertEqual(g.AtomCharges, lsg.AtomCharges, "layersubgeometry charges differ")
        
    
    def testXYZio(self):
        """.xyz i/o"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.xyz")
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
        g2=comatsci.geometry.Geometry()
        g2.readfile("test.xyz")
        os.unlink("test.xyz")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        # append mode
        self.assertRaises(comatsci.geometry.GeometryError, g.writexyz, "test.xyz", mode="x") # invalid write mode
        g.writexyz("test.xyz")
        g.writexyz("test.xyz",mode="a")
        g2=comatsci.geometry.Geometry()
        self.assertRaises(ValueError, g2.readfile, "test.xyz") # multi frame xyz not supported
        testfile=open("test.xyz")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.xyz")
        g2.parseXyzString("".join(testlines[len(testlines)/2:]))
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
    
    
    def testGENio(self):
        """.gen i/o"""
        # cluster geometries
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.gen")
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
        g2=comatsci.geometry.Geometry()
        g2.readfile("test.gen")
        os.unlink("test.gen")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        self.assertRaises(comatsci.geometry.GeometryError, g.writegen, "test.gen", cmode="B") # invalid coordinate mode
        self.assertRaises(comatsci.geometry.GeometryError, g.writegen, "test.gen", cmode="F") # fractional coordinates illegal for cluster geometry
        # periodic
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest2.gen")
        g2=comatsci.geometry.Geometry()
        self.assertEqual(g.Mode, "S", "periodic geometry mode differs")
        georef=numpy.array([[0.000000E+00,             0.000000E+00,             0.000000E+00 ],
          [  4.340446E+00 ,            4.340446E+00 ,            2.795018E+00 ],
          [  2.650276E+00 ,            2.650276E+00 ,            0.000000E+00 ],
          [ -2.650276E+00 ,           -2.650276E+00 ,            0.000000E+00 ],
          [  6.990722E+00 ,            1.690170E+00 ,            2.795018E+00 ],
          [  1.690170E+00 ,            6.990722E+00 ,            2.795018E+00]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry coordinates differ")
        latticeref=numpy.array([[ 8.68089177,  0.        ,  0.        ],  # @UnusedVariable
                          [ 0.        ,  8.68089177,  0.        ],
                          [ 0.        ,  0.        ,  5.59003676]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry lattice vectors differ")
        g.writegen("test.gen")
        g2.readgen("test.gen")
        os.unlink("test.gen")
        self.assertEqual("S",g2.Mode,"written periodic geometry mode differs")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written periodic geometry coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "written periodic geometry lattice vectors differ")
        # fractional coordinates
        g.writegen("test.gen",cmode="F")
        g2.readgen("test.gen")
        os.unlink("test.gen")
        self.assertEqual("S",g2.Mode,"written fractional coordinates geometry mode differs")
        print(g.Geometry,g2.Geometry)
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written fractional coordinates geometry coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "written fractional coordinates geometry lattice vectors differ")
    
    
    
    def testCDHio(self):
        """.cdh i/o"""
        g=comatsci.geometry.Geometry()
        g.readCDHFile("../../test_regr/unittest_data/geotest1.cdh")
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
        g2=comatsci.geometry.Geometry()
        g2.readCDHFile("test.cdh")
        os.unlink("test.cdh")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        # periodic
        g=comatsci.geometry.Geometry()
        g.readCDHFile("../../test_regr/unittest_data/geotest2.cdh")
        g2=comatsci.geometry.Geometry()
        self.assertEqual(g.Mode, "S", "periodic geometry mode differs")
        georef=numpy.array([[0.000000E+00,             0.000000E+00,             0.000000E+00 ],
          [  4.340446E+00 ,            4.340446E+00 ,            2.795018E+00 ],
          [  2.650276E+00 ,            2.650276E+00 ,            0.000000E+00 ],
          [ -2.650276E+00 ,           -2.650276E+00 ,            0.000000E+00 ],
          [  6.990722E+00 ,            1.690170E+00 ,            2.795018E+00 ],
          [  1.690170E+00 ,            6.990722E+00 ,            2.795018E+00]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry coordinates differ")
        latticeref=numpy.array([[ 8.68089177,  0.        ,  0.        ],  # @UnusedVariable
                          [ 0.        ,  8.68089177,  0.        ],
                          [ 0.        ,  0.        ,  5.59003676]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry lattice vectors differ")
        g.writeCDH("test.cdh")
        g2.readCDHFile("test.cdh")
        os.unlink("test.cdh")
        self.assertEqual("S",g2.Mode,"written periodic geometry mode differs")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written periodic geometry coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "written periodic geometry lattice vectors differ")
        # globals
        g2=comatsci.geometry.Geometry()
        g2.readCDHFile("../../test_regr/unittest_data/geotest2b.cdh")
        self.assertEqual(g.Mode,g2.Mode,"globals CDH file mode differs")
        self.assertEqual(g.AtomTypes,g2.AtomTypes,"globals CDH file elements differ")
        self.assertEqual(g.AtomSubTypes,g2.AtomSubTypes,"globals CDH file types differ")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "globals CDH file coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "globals CDH file lattice vectors differ")
        
    
    
    
    def testFMGio(self):
        """.fmg i/o"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
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
        g2=comatsci.geometry.Geometry()
        g2.readfile("test.fmg")
        os.unlink("test.fmg")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written geometry coordinates differ")
        self.assertEqual([6,6,1,1,1,1], g2.AtomTypes, "written geometry types differ")
        self.assertEqual(["C","C","H","H","H","H"], g2.AtomSubTypes, "written geometry subtypes differ")
        # periodic
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest2.fmg")
        g2=comatsci.geometry.Geometry()
        self.assertEqual(g.Mode, "S", "periodic geometry mode differs")
        georef=numpy.array([[0.000000E+00,             0.000000E+00,             0.000000E+00 ],
          [  4.340446E+00 ,            4.340446E+00 ,            2.795018E+00 ],
          [  2.650276E+00 ,            2.650276E+00 ,            0.000000E+00 ],
          [ -2.650276E+00 ,           -2.650276E+00 ,            0.000000E+00 ],
          [  6.990722E+00 ,            1.690170E+00 ,            2.795018E+00 ],
          [  1.690170E+00 ,            6.990722E+00 ,            2.795018E+00]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry coordinates differ")
        latticeref=numpy.array([[ 8.68089177,  0.        ,  0.        ],  # @UnusedVariable
                          [ 0.        ,  8.68089177,  0.        ],
                          [ 0.        ,  0.        ,  5.59003676]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry lattice vectors differ")
        g.writegen("test.fmg")
        g2.readgen("test.fmg")
        os.unlink("test.fmg")
        self.assertEqual("S",g2.Mode,"written periodic geometry mode differs")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written periodic geometry coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "written periodic geometry lattice vectors differ")
        
        
    def testAimsOut(self):
        """fhi-aims i/o"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
        g.writeAIMS("test.aims")
        testfile=open("test.aims")
        testlines=list(testfile)
        testfile.close()
        #read
        g2=comatsci.geometry.Geometry()
        g2.readAimsFile("test.aims")
        os.unlink("test.aims")
        self.assertEqual(g.getAimsString(), "".join(testlines).strip("\n"), "written file corrupted")
        AtomTypes=["C","C","H","H","H","H"]
        georef=numpy.array([[ 0.          ,0.      ,    0.        ],
                            [ 0.          ,0.      ,    2.52278443],
                            [ 1.78220452  ,0.      ,   -1.0289559 ],
                            [-1.78220452  ,0.      ,    3.55174032],
                            [ 1.78220452  ,0.      ,    3.55174032],
                            [-1.78220452  ,0.      ,   -1.0289559 ]],dtype=float)
        georef*=comatsci.geometry.Angstrom
        for ii in range(g.Atomcount):
            tokens=testlines[ii+1].split()
            self.assertEqual(AtomTypes[ii],tokens[4],"Element type differs")
            self.assertAlmostEqual(georef[ii][0],float(tokens[1]),6,"X coordinate differs")
            self.assertAlmostEqual(georef[ii][1],float(tokens[2]),6,"Y coordinate differs")
            self.assertAlmostEqual(georef[ii][2],float(tokens[3]),6,"Z coordinate differs")
        # test read
        self.assertEqual(g.Mode,g2.Mode,"read aims geometry mode differs")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "read aims geometry coordinates differ")
        self.assertEqual(AtomTypes,g2.AtomSubTypes,"read aims atom types differ")
        # periodic
        g=comatsci.geometry.Geometry()
        g.readAimsFile("../../test_regr/unittest_data/geotest2.aims")
        g2=comatsci.geometry.Geometry()
        self.assertEqual(g.Mode, "S", "periodic geometry mode differs")
        georef=numpy.array([[0.000000E+00,             0.000000E+00,             0.000000E+00 ],
          [  4.340446E+00 ,            4.340446E+00 ,            2.795018E+00 ],
          [  2.650276E+00 ,            2.650276E+00 ,            0.000000E+00 ],
          [ -2.650276E+00 ,           -2.650276E+00 ,            0.000000E+00 ],
          [  6.990722E+00 ,            1.690170E+00 ,            2.795018E+00 ],
          [  1.690170E+00 ,            6.990722E+00 ,            2.795018E+00]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry coordinates differ")
        latticeref=numpy.array([[ 8.68089177,  0.        ,  0.        ],  # @UnusedVariable
                          [ 0.        ,  8.68089177,  0.        ],
                          [ 0.        ,  0.        ,  5.59003676]],dtype=float)
        self.assertTrue(numpy.allclose(georef, g.Geometry), "periodic geometry lattice vectors differ")
        g.writegen("test.aims")
        g2.readgen("test.aims")
        os.unlink("test.aims")
        self.assertEqual("S",g2.Mode,"written periodic geometry mode differs")
        self.assertTrue(numpy.allclose(g.Geometry, g2.Geometry), "written periodic geometry coordinates differ")
        self.assertTrue(numpy.allclose(g.Lattice, g2.Lattice), "written periodic geometry lattice vectors differ")
        


    def testTmOut(self):
        """turbomole output"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
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
        # periodic
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest2.fmg")
        self.assertRaises(comatsci.geometry.GeometryError, g.writeTurboMole, "test.tm")
        
    
    
    def testFdfOut(self):
        """SIESTA .fdf output"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
        g.writefdf("test.fdf")
        testfile=open("test.fdf")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.fdf")
        reffile=open("../../test_regr/unittest_data/geotest1.fdf")
        reflines=list(reffile)
        reffile.close()
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written fdf file differs from reference")
        #periodic
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest2.fmg")
        g.writefdf("test.fdf")
        testfile=open("test.fdf")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.fdf")
        reffile=open("../../test_regr/unittest_data/geotest2.fdf")
        reflines=list(reffile)
        reffile.close()
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written periodic fdf file differs from reference")
        
        


    def testPDBOut(self):
        """.pdb output"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
        g.writepdb("test.pdb")
        testfile=open("test.pdb")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.pdb")
        reffile=open("../../test_regr/unittest_data/geotest1.pdb")
        reflines=list(reffile)
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written pdb file differs from reference")
        #periodic
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest2.fmg")
        g.writepdb("test.pdb")
        testfile=open("test.pdb")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.pdb")
        reffile=open("../../test_regr/unittest_data/geotest2.pdb")
        reflines=list(reffile)
        reffile.close()
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written periodic .pdb file differs from reference")
        #errors
        self.assertRaises(comatsci.geometry.GeometryError, g.writepdb, "test.pdb", beta=[1,2,3])
        self.assertRaises(comatsci.geometry.GeometryError, g.writepdb, "test.pdb", occupancy=[1,2,3])
        # bond list
        g.writepdb("test.pdb",writebondlist=1)
        testfile=open("test.pdb")
        testlines=list(testfile)
        testfile.close()
        os.unlink("test.pdb")
        reffile=open("../../test_regr/unittest_data/geotest2b.pdb")
        reflines=list(reffile)
        reffile.close()
        for i in range(len(testlines)):
            self.assertEqual(reflines[i],testlines[i],"written periodic connect records .pdb file differs from reference")
        
            
            
    def testXYZQOut(self):
        """DFTB+ .xyzq output"""
        g=comatsci.geometry.Geometry()
        g.readfile("../../test_regr/unittest_data/geotest1.fmg")
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
    