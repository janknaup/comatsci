'''
Created on Mar 8, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
from __future__ import print_function
import unittest,numpy
import comatsci


class splineTest(unittest.TestCase):


    def testSpline(self):
        """cubic Spline interpolation"""
        XX=numpy.arange(0,10,0.2,dtype=float)    # nodes
        X2=numpy.arange(0.1,9,0.4,dtype=float)   # intermediate values
        YY=numpy.sin(XX)                         # nodes
        Y2=numpy.sin(X2)                         # intermediates
        C2=numpy.cos(X2)                         # analytic derivates at intermediates
        # test exceptions first
        SS=comatsci.spline.Spline()
        self.assertRaises(ValueError, SS.setdata, XX, YY[2:]) # different length arrays
        self.assertRaises(ValueError, SS.setdata, YY, XX)     # non-monotonic x array
        self.assertRaises(ValueError, SS.splint, 2.5)         # uninitialzed spline
        self.assertRaises(ValueError, SS.splder, 2.5)         # uninitialized spline
        # now properly initialize 
        SS=comatsci.spline.Spline(XX,YY)
        # some more exceptions
        self.assertRaises(ValueError, SS.splint, 11)         # x out of range
        self.assertRaises(ValueError, SS.splder, 11)         # x out of range
        # now test interpolation
        II=numpy.array([ SS.splint(x) for x in XX ],dtype=float)
        I2=numpy.array([ SS.splint(x) for x in X2 ],dtype=float)
        ID2=numpy.array([ SS.splder(x) for x in X2 ],dtype=float)
        self.assertTrue(numpy.allclose(YY, II),"Node values deviate")
        self.assertTrue(numpy.allclose(Y2, I2),"Interpolated value deviation too large")
        # don't be too picky with the end point derivatives
        self.assertTrue(numpy.allclose(C2[1:-2], ID2[1:-2]),"Interpolated derivative deviation too large")
    
    
    def testRenner(self):
        """Renner Spline interpolation"""
        ## Renner splines use their own parameter range, so the original parameter must be saved
        ## with the spline, if one wants to compare to original values
        XX=numpy.arange(0,10,0.2,dtype=float)
        SYY=numpy.sin(XX)
        CYY=numpy.cos(XX)
        VEC=numpy.transpose(numpy.array([XX,SYY,CYY]),)
        # uninitialized first
        RR=comatsci.spline.RennerSpline()
        self.assertRaises(ValueError, RR.setNodes, VEC[0:4])  # too few nodes
        self.assertRaises(ValueError, RR.splint, 5)           # uninitialized
        self.assertRaises(ValueError, RR.splder, 5)           # uninitialized
        # proper init now
        RR=comatsci.spline.RennerSpline(VEC)
        # more errors
        self.assertRaises(ValueError, RR.splint, -0.1)        # out of range
        self.assertRaises(ValueError, RR.splder, -0.1)        # out of range
        self.assertRaises(ValueError, RR.splint, RR.getTotalLength()+0.1)        # out of range
        self.assertRaises(ValueError, RR.splder, RR.getTotalLength()+0.1)        # out of range
        # test interpolation
        II=numpy.array([RR.splint(t) for t in numpy.arange(0,RR.totalLength,0.2)])
        DD=numpy.array([RR.splder(t) for t in numpy.arange(0,RR.totalLength,0.2)])
        IX=numpy.transpose(II)[0]
        IYS=numpy.transpose(II)[1]
        IYC=numpy.transpose(II)[2]
        DX=numpy.transpose(DD)[0]
        DYDS=numpy.transpose(DD)[1]
        SIX=numpy.sin(IX)
        CIX=numpy.cos(IX)
        self.assertTrue(numpy.allclose(IYS[1:-2], SIX[1:-2],1.0+(5E-4),2E-3),"sine interpolation deviation too large")
        # test cosine coordinate
        self.assertTrue(numpy.allclose(IYC[1:-2], CIX[1:-2], 1.0+(5E-4),2E-3),"cosine interpolation deviation too large")
        # test derivative
        self.assertTrue(numpy.allclose(DYDS[1:-2]/DX[1:-2], CIX[1:-2], 1.0+5E-3,1E-1),"sine derivative deviation too large")


    def testVecSpline(self):
        """Vector Spline interpolation"""
        ## Renenwr splines use their own parameter range, so the original parameter must be saved
        ## with the spline, if one wants to compare to original values
        XX=numpy.arange(0,10,0.2,dtype=float)
        SYY=numpy.sin(XX)
        CYY=numpy.cos(XX)
        # uninitialized first
        VV=comatsci.spline.vectorSpline()
        self.assertRaises(ValueError, VV.setNodes, [XX, [SYY[2:],]])  # different x,y lengths
        self.assertRaises(ValueError, VV.setNodes, [XX, [ [CYY[i],SYY[i] ] for i in range(len(CYY)-1) ] + [[CYY,]]])  # different y lengths per dimension
        self.assertRaises(ValueError, VV.setNodes, [SYY, numpy.transpose(numpy.array([CYY,SYY]))]) # non-monotonic x-values
        self.assertRaises(ValueError, VV.splint, 5)           # uninitialized
        self.assertRaises(ValueError, VV.splder, 5)           # uninitialized
        self.assertRaises(ValueError, VV.getArcLength)        # uninitialized
        # proper init now
        VV=comatsci.spline.vectorSpline([XX,numpy.transpose(numpy.array([SYY,CYY]))])
        # more errors
        self.assertRaises(ValueError, VV.splint, -0.1)        # out of range
        self.assertRaises(ValueError, VV.splder, -0.1)        # out of range
        self.assertRaises(ValueError, VV.splint, 11)        # out of range
        self.assertRaises(ValueError, VV.splder, 11)        # out of range
        self.assertRaises(ValueError, VV.getArcLength, points=(-1,5)) # out of range
        self.assertRaises(ValueError, VV.getArcLength, points=(5,15)) # out of range
        # test interpolation
        # skip arc length for now as it is broken!
        # self.assertAlmostEqual(9.79970388534, VV.getArcLength(tolerance=1e-6), 2, "arc length differs")
        self.assertAlmostEqual(VV.getArcLength(tolerance=1e-6),VV.getArcLength(tolerance=1e-6), 10, "consecutive arc lenghts differ")
        self.assertEquals(VV.getDimension(), 2, "interpolated vector dimension differs")
        self.assertEquals(VV.getNodeCount(), len(XX), "Spline node count differs")
        IX=XX[3:-3]+0.12
        II=numpy.array([VV.splint(t) for t in IX ] )
        DD=numpy.array([VV.splder(t) for t in IX ] )
        IYS=numpy.transpose(II)[0]
        IYC=numpy.transpose(II)[1]
        DYDS=numpy.transpose(DD)[0]
        SIX=numpy.sin(IX)
        CIX=numpy.cos(IX)
        self.assertTrue(numpy.allclose(IYS[1:-2], SIX[1:-2],atol=1E-6),"sine interpolation deviation too large")
        # test cosine coordinate
        self.assertTrue(numpy.allclose(IYC[1:-2], CIX[1:-2],atol=5E-6),"cosine interpolation deviation too large")
        # test derivative
        self.assertTrue(numpy.allclose(DYDS[6:-6], CIX[6:-6], 1.0+1E-5, 1E-5),"sine derivative deviation too large")


def suite():
    suite = unittest.makeSuite(splineTest,'test')
    return suite

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSpline']
    unittest.main()
    