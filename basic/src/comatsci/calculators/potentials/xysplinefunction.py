##############################################################################
# xySplineFunction.py
# Part of comatsci computational materials science toolkit
# (c) 2011 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
# adapted from polynomFunction by Manuel Doemer <manuel.doemer@epfl.ch>
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <janknaup@gmail.com>
#@organization: EPFL-ISIC-LCBC
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup  <janknaup@gmail.com>

from __future__ import print_function
from .potentialfunction import potentialFunction
from ..calcerror import CalcError
from ... import Spline

import numpy 

class xySplineFunction(potentialFunction):
		""" class for r,E(r) Spline represenation of a pairpotential """

		def __init__(self,parameters):
				""" construct r,E(r) Spline
				@parameter parameters: dictionary of function parameters """

				potentialFunction.__init__(self,parameters)

				if self._parameters.has_key("x") and self._parameters.has_key("y"):
					if len(parameters["x"])!=len(parameters["y"]):
						raise ValueError("x and y arrays of xy-spline potential not equal")
					self._spline=Spline(self._parameters.get("x"),self._parameters.get("y"))



		def parseString(self,inString):
				""" parse spline x,y values from inString string
				update self._parameters
				@parameter inString: string containing function parameters, separated by \n """
				inlines=inString.split("\n")
				""" first line should contain &xyspline& """
				dummy=inlines[0].split()[0].strip()
				if not (dummy=='&xyspline&'):
						raise CalcError(" format of xy spline file not recognized ")
				# remember: the first line should contain the potential function identifier field, e.g. &polynom&
				xvals=[]
				yvals=[]
				boundaries=[]
				mutables=[]
				for i in range(1,len(inlines)):
					# ignore common comment markers and blank lines
					if len(inlines[i])==0:
						firstchar="#"
					else:
						firstchar=inlines[i].lstrip()[0]
					if not (firstchar=="#" or firstchar==";" or firstchar=="!" or firstchar=='*'):
						dummy=inlines[i].split()
						try:
							xvals.append(float(dummy[0]))
							yvals.append(float(dummy[1]))
						except:
							print("***ERROR: Error parsing xy spline inString file line {0:d}".format(i+1))
							raise
						if len(dummy)==2:
							boundaries.append(None)
						elif len(dummy)==3:
							if dummy[2].lower() in ("true","t","y","1","yes"):
								mutables.append(len(xvals)-1)
							boundaries.append(None)
						elif len(dummy)==4:
							try:
								boundaries.append(
									(
										float(dummy[2]),
										float(dummy[3])
									)
								)
							except:
								print("***ERROR: Error parsing xy spline inString file line {0:d}".format(i+1))
								raise
							mutables.append(len(xvals)-1)
						else:
							print("***WARNING: trailing garbage in xy spline inString file")
				# store inner and outer cutoff values
				self._innercut=xvals[0]
				self._outercut=xvals[-1]
				self._boundaries=boundaries
				self._mutables=mutables
				#store parameters
				self._parameters={
					"x": numpy.array(xvals),
					"y": numpy.array(yvals)
				}
				# reinit spline rep
				self._spline=Spline(self._parameters["x"],self._parameters["y"])
				#done.

		
		def value(self,r):
				""" calculates the value of the spline
				at a given distance r 
				@parameter r: float containing distance 
				@return: float containing function value """
				value=None
				if (r < self._innercut) or (r > self._outercut):
					raise ValueError("Interpolation requested for out-of-range x value")
				else:
					value=self._spline.splint(r)
				return value


		def derivative(self,r):
				""" returns derivative of the polynom
				at a given distance r
				@parameter r: float containing distance 
				@return: float containing function derivative value """
				
				derivative=None
				if (r < self._innercut) or (r > self._outercut):
					raise ValueError("Interpolation requested for out-of-range x value")
				else:
					derivative=self._spline.splder(r)
				return derivative


		def getCutoffs(self):
				""" @return: tuple with inner and outer cutoff """
				return (self._innercut,self._outercut)
		
		
		
		
		def getMutables(self):
				""" @return: list with all parameters allowed to change """
				mutlist=[]
				for i in self._mutables:
					mutlist.append(self._parameters["y"][i])
				return mutlist
			
			
		
		def updateMutables(self,mutated):
				""" @parameter mutated: list with updated parameters to write back
				into self._parameters
				"""
				# check if number of passed values matches mutables
				if len(mutated)!=len(self._mutables):
					raise ValueError("Number of mutated values does not match number of mutables in potential")
				# copy mutated values into parameters dictionary
				for i in range(len(mutated)):
					self._parameters["y"][self._mutables[i]]=mutated[i]
				# update internal spline representation
				self._spline=Spline(self._parameters["x"],self._parameters["y"])
				#done.




		def getString(self):
			"""@return: string containing potential function file representation of the current instance"""
			potLines=["&xyspline&"]
			for i in range(len(self._parameters["x"])):
				line="{0:12f}\t{1:12f}".format(self._parameters["x"][i],self._parameters["y"][i])
				if i in self._mutables:
					if self._boundaries[i]!=(None,None):
						line+="\t{0[0]:12f}\t{0[1]:12f}".format(self._boundaries[i])
					else:
						line+="\tTrue"
				potLines.append(line)
			return "\n".join(potLines)
		


