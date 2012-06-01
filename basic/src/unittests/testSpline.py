'''
Created on Mar 8, 2012

@author: Jan M. Knaup <janknaup@gmail.com>
'''
import unittest,numpy
from comatsci import Spline


class splineTest(unittest.TestCase):


    def testSpline(self):
        """Test cubic Spline interpolation"""
        XX=numpy.arange(0,10,0.2,dtype=float)    # nodes
        X2=numpy.arange(0.1,9,0.4,dtype=float)   # intermediate values
        YY=numpy.sin(XX)                         # nodes
        Y2=numpy.sin(X2)                         # intermediates
        C2=numpy.cos(X2)                         # analytic derivates at intermediates
        SS=Spline.spline(XX,YY)
        II=numpy.array([ SS.splint(x) for x in XX ],dtype=float)
        I2=numpy.array([ SS.splint(x) for x in X2 ],dtype=float)
        ID2=numpy.array([ SS.splder(x) for x in X2 ],dtype=float)
        self.assertTrue(numpy.allclose(YY, II),"Node values deviate")
        self.assertTrue(numpy.allclose(Y2, I2),"Interpolated value deviation too large")
        # don't be too picky with the end point derivatives
        self.assertTrue(numpy.allclose(C2[1:-2], ID2[1:-2]),"Interpolated derivative deviation too large")
    
    
    def testRenner(self):
        """Test Renner Spline interpolation"""
        ## Renenwr splines use their own parameter range, so the original parameter must be saved
        ## with the spline, if one wants to compare to original values
        XX=numpy.arange(0,10,0.2,dtype=float)
        SYY=numpy.sin(XX)
        CYY=numpy.cos(XX)
        VEC=numpy.transpose(numpy.array([XX,SYY,CYY]),)
        RR=Spline.RennerSpline(VEC)
        II=numpy.array([RR.splint(t) for t in numpy.arange(0,RR.totalLength,0.2)])
        IX=numpy.transpose(II)[0]
        IYS=numpy.transpose(II)[1]
        IYC=numpy.transpose(II)[2]
        SIX=numpy.sin(IX)
        CIX=numpy.cos(IX)
        self.assertTrue(numpy.allclose(IYS[1:-2], SIX[1:-2],5E-4,2E-3),"sine interpolation deviation too large")
        # test cosine coordinate
        self.assertTrue(numpy.allclose(IYC[1:-2], CIX[1:-2], 5E-4,2E-3),"cosine interpolation deviation too large")
#        pf=open("renner.dat","w")
#        for ii in range(len(IX)):
#            print >> pf, IX[ii], IYS[ii], SIX[ii], IYC[ii], CIX[ii]
#        pf.close() 


def suite():
    suite = unittest.makeSuite(splineTest,'test')
    return suite

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSpline']
    unittest.main()