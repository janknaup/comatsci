##############################################################################
# potentialFunction.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <Knaup@bccms.uni-bremen.de>
#@organization: Bremen Center for Compuational Materials Science
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup  <Knaup@bccms.uni-bremen.de>

from comatsci.Calculators.CalcError import *

class potentialFunction:
        """ base class for pairpotential functions """

        def __init__(self,parameters=None):
                """ pairpotential constructor
                @parameter parameters: dictionary of funtion parameters """
                if (parameters==None):
                        self._parameters=dict()
                elif not type(parameters) is dict:
                        raise CalcError("parameters is not of type dictionary")
                else:        
                        self._parameters=parameters
                # the boundaries dictionary sets boundaries to toe parameters
                # being subject to optimization
                self._boundaries=self._parameters.pop("boundaries",{})



        def value(self,r):
                """ calculates the value of the respective pairpotential
                at a given distance r 
                @parameter r: float containing distance 
                @return float containing function value """

                raise NotImplementedError()



        def derivative(self,r):
                """ returns derivative of the respective pairpotential
                at a given distance r
                @parameter r: float containing distance 
                @return float containing function derivative value """

                raise NotImplementedError()



        def getParameters(self):
                """ @return dictionary of function parameters """
                
                return dict(self._parameters)



        def readFile(self,filename):
                """ read potential parameters from file 
                @parameter filename: string containing input file name"""

                infile=open(filename,"r")
                inlines=list(infile)
                infile.close()
                self.parseString("".join(inlines))



        def writeFile(self,filename):
                """ write potential parameters to file
                @parameter filename: string containing output file name"""
                
                outfile=open(filename,"w")
                outfile.write(self.getString())
                outfile.close()



        def parseString(self,input):
                """ parse string containing function parameters
                @parameter input: string containing function parameters """

                raise  NotImplementedError()



        def getString(self):
                """ convert parameter object to string """

                raise  NotImplementedError()


        def getCutoffs(self):
                """ @return cutoffs as a tuple """

                raise NotImplementedError()

        def updateMutables(self):
                """ @return list with all parameters allowed to change """
                raise NotImplementedError()

        def updateMutables(self,mutated):
                """ @parameter mutated: list with updated parameters to write back
                into self._parameters
                """

                raise NotImplementedError()




