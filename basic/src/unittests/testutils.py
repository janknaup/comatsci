'''
Created on Oct 30, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
from __future__ import print_function
import unittest,numpy,os
import comatsci

class UtilsTest(unittest.TestCase):
    """test comatsci utilities"""
    
    
    def testCompressOpen(self):
        """compressed file opening"""
        self.assertRaises(ValueError, comatsci.utils.compressedopen, "doesnotexist", mode="r", compresslevel=0, autodetect=True)
        self.assertRaises(ValueError, comatsci.utils.compressedopen, "doesnotexist", mode="r", compresslevel=0, autodetect=False)
        self.assertRaises(ValueError, comatsci.utils.compressedopen, "doesnotexist", mode="r", compresslevel=-2, autodetect=False)
        self.assertRaises(ValueError, comatsci.utils.compressedopen, "doesnotexist", mode="r", compresslevel=10, autodetect=False)
        reffile=open("../../test_regr/unittest_data/geotest1.gen")
        reflines=list(reffile)
        reffile.close()
        gzfile=comatsci.utils.compressedopen("../../test_regr/unittest_data/geotest1b.gen.gz")
        gzlines=list(gzfile)
        gzfile.close
        for ii in range(len(reflines)):
            self.assertEqual(reflines[ii], gzlines[ii], "read gzip lines differ")
        bzfile=comatsci.utils.compressedopen("../../test_regr/unittest_data/geotest1b.gen.bz2")
        bzlines=list(bzfile)
        bzfile.close
        for ii in range(len(reflines)):
            self.assertEqual(reflines[ii], bzlines[ii], "read bzip2 lines differ")
        self.assertRaises(ValueError, comatsci.utils.compressedopen, "../../test_regr/unittest_data/geotest1b.gen", autodetect=False)
        bzfile=comatsci.utils.compressedopen("../../test_regr/unittest_data/geotest1b.gen", autodetect=True)
        bzlines=list(bzfile)
        bzfile.close()
        gzfile=comatsci.utils.compressedopen("test.gen.gz",mode="w",compresslevel=9)
        for ii in range(len(reflines)):
            gzfile.write(reflines[ii])
        gzfile.close()
        gzfile=comatsci.utils.compressedopen("test.gen", autodetect=True)
        gzlines=list(gzfile)
        gzfile.close()
        for ii in range(len(reflines)):
            self.assertEqual(reflines[ii], gzlines[ii], "written gzip lines differ")
        os.unlink("test.gen.gz")
        # writing bz2 file
        gzfile=comatsci.utils.compressedopen("test.gen.bz2",mode="w",compresslevel=9)
        for ii in range(len(reflines)):
            gzfile.write(reflines[ii])
        gzfile.close()
        gzfile=comatsci.utils.compressedopen("test.gen", autodetect=True)
        gzlines=list(gzfile)
        gzfile.close()
        for ii in range(len(reflines)):
            self.assertEqual(reflines[ii], gzlines[ii], "written bzip2 lines differ")
        os.unlink("test.gen.bz2")
        # autodetected files for writing
        gzfile=comatsci.utils.compressedopen("test.dat", mode="w", compresslevel=0)
        print("test data",file=gzfile)
        gzfile.close()
        self.assertTrue(os.path.exists("test.dat"), "uncompressed file was not created")
        gzfile=comatsci.utils.compressedopen("test.dat","r")
        self.assertEqual(gzfile.readline(), "test data\n", "written uncompressed file contents differ")
        gzfile.close()
        os.unlink("test.dat")
        gzfile=comatsci.utils.compressedopen("test.dat", mode="w", compresslevel=5)
        print("test data",file=gzfile)
        gzfile.close()
        self.assertTrue(os.path.exists("test.dat.gz"), "compressed file was not created")
        gzfile=comatsci.utils.compressedopen("test.dat","r")
        self.assertEqual(gzfile.readline(), "test data\n", "written uncompressed file contents differ")
        gzfile.close()
        os.unlink("test.dat.gz")
        
    
    def testUncompressCopy(self):
        """copy and uncompress"""
        # gzipped
        testfile=comatsci.utils.compressedopen("test.txt.gz", mode="w", compresslevel=9)
        print("test text",file=testfile,end="")
        testfile.close()
        comatsci.utils.uncompresscopy("test.txt.gz", "copied.txt")
        testfile=open("copied.txt")
        testlines=list(testfile)
        testfile.close()
        self.assertEqual("test text", testlines[0], "uncompressed file contents differ")
        os.unlink("copied.txt")
        os.unlink("test.txt")
        # bzip2
        testfile=comatsci.utils.compressedopen("test.txt.bz2", mode="w", compresslevel=9)
        print("test text",file=testfile,end="")
        testfile.close()
        comatsci.utils.uncompresscopy("test.txt.bz2", "copied.txt")
        testfile=open("copied.txt")
        testlines=list(testfile)
        testfile.close()
        self.assertEqual("test text", testlines[0], "uncompressed file contents differ")
        os.unlink("copied.txt")
        os.unlink("test.txt")
        # plain
        testfile=comatsci.utils.compressedopen("test.txt", mode="w", compresslevel=0)
        print("test text",file=testfile,end="")
        testfile.close()
        comatsci.utils.uncompresscopy("test.txt", "copied.txt")
        testfile=open("copied.txt")
        testlines=list(testfile)
        testfile.close()
        self.assertEqual("test text", testlines[0], "copied file contents differ")
        os.unlink("copied.txt")
        os.unlink("test.txt")
        
        
    
    
    def testCompressCopy(self):
        """compress and copy"""
        # copy compressed
        testfile=open("test.txt","w")
        print("test text",file=testfile,end="")
        testfile.close()
        comatsci.utils.compresscopy("test.txt", "copied.txt.gz")
        #self.assertTrue(os.path.exists("test.txt"), "original file missing")
        self.assertTrue(os.path.exists("copied.txt.gz"), "copied compressed file missing")
        testfile=comatsci.utils.compressedopen("copied.txt.gz")
        testlines=list(testfile)
        testfile.close()
        self.assertEqual("test text", testlines[0], "copied compressed file contents differ")
        self.assertRaises(ValueError, comatsci.utils.compresscopy, "test.txt", "copied.txt.gz", compresslevel=0)
        self.assertRaises(ValueError, comatsci.utils.compresscopy, "test.txt", "copied.txt.gz", compresslevel=10)
        os.unlink("copied.txt.gz")
    
    
    def testProgressMeter(self):
        """progress meter"""
        # just test for Errors
        pm=comatsci.utils.ProgressMeter()
        for i in range(0,101):
            pm.update(i)
    
    
    def testPrettyPrint(self):
        """dictionary pretty printer"""
        dict1={"eins":1,"zwei":2,"drei":3}
        dict2={"zahlen":dict1,"text":"text","float":3.14157}
        string1=comatsci.utils.dictionaryPrettyPrint(dict1)
        string2=comatsci.utils.dictionaryPrettyPrint(dict2)
        self.assertEqual(string1,"drei: 3\neins: 1\nzwei: 2","single level dictionary string differs")
        self.assertEqual(string2,"float : 3.14157\ntext  : text\nzahlen: \n         drei: 3\n         eins: 1\n         zwei: 2","two level dictionary string differs")
    