##############################################################################
# slaterFunction.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Manuel Doemer <manuel.doemer@epfl.ch>
#@organization: Ecole Polytechnique Federale de Lausanne
#@license: Open Software License version 3.0
#@copyright: Manuel Doemer	<manuel.doemer@epfl.ch>

from .potentialfunction import potentialFunction
from ..calcerror import CalcError
import numpy
import math


class slaterFunction(potentialFunction):
	""" class for slater represenation of a pairpotential """

	def __init__(self, parameters):
		""" slater constructor
		@parameter parameters: dictionary of function parameters """

		potentialFunction.__init__(self, parameters)



	def parseString(self, inString):
		""" read slater coefficients from inString string
		update self._parameters
		@parameter inString: string containing function parameters, separated by \n """
		inlines = inString.split("\n")
		""" first line should contain &slater& """
		dummy = inlines[0].split()[0].strip()
		if not (dummy == '&slater&'):
			raise CalcError(" format of slater coefficient file not recognized ")

		dummy = ''
		counter = 0 # counts the number of lines read-in as containing parameters

		for i in inlines[1:len(inlines)]:
			# ignore common comment markers and blank lines 
			if not (i == ''):
				firstchar = i.lstrip()[0]
				if not (firstchar == "#" or firstchar == ";" or firstchar == "!" or firstchar == '*'):
					dummy = i.split()
					counter = counter + 1
					# this counts also the &potentialFunction& line !
					# the first parameters containing line holds the cutoff
					if (counter == 1):
						self._parameters['cutoff'] = float(dummy[0])
						if (len(dummy) == 3):
							# if there are two more numbers in the line
							# the cutoff is assumed to be subject to optimization
							# within the given boundaries
							self._boundaries['cutoff'] = (float(dummy[1]), float(dummy[2]))
						elif (len(dummy) == 1):
							# no boundaries given: cutoff will be held fixed
							# during optimization
							self._boundaries['cutoff'] = None
						else:
							raise CalcError(" incorrect no of variable cutoff boundaries ")
					else:
						if (int(dummy[0]) < 1):
							raise CalcError(" parameter n < 1 !")

						else:
							self._parameters['n' + str(counter - 1)] = int(dummy[0])
						
						if (len(dummy) == 3):
							# c_i and z_i are read, but assumed to be kept fixed
							self._parameters['c' + str(counter - 1)] = float(dummy[1])
							self._boundaries['c' + str(counter - 1)] = None
							self._parameters['z' + str(counter - 1)] = float(dummy[2])
							self._boundaries['z' + str(counter - 1)] = None
						
						elif (len(dummy) == 5):
							# c_i is read and assumed to be subject to optimzation
							# within the given boundaries
							# z_i is read as fifth field and assumed to be kept fixed
							self._parameters['c' + str(counter - 1)] = float(dummy[1])
							self._boundaries['c' + str(counter - 1)] = (float(dummy[2]), float(dummy[3]))
							self._parameters['z' + str(counter - 1)] = float(dummy[4])
							self._boundaries['z' + str(counter - 1)] = None

						elif (len(dummy) == 7):
							# c_i and z_i are read with the respective boundaries
							# and are assumed to be subject to optimization
							self._parameters['c' + str(counter - 1)] = float(dummy[1])
							self._boundaries['c' + str(counter - 1)] = (float(dummy[2]), float(dummy[3]))
							self._parameters['z' + str(counter - 1)] = float(dummy[4])
							self._boundaries['z' + str(counter - 1)] = (float(dummy[5]), float(dummy[6]))

						else:
							raise CalcError(" incorrect no of variable parameter boundaries ")
		# number of slater functions in the expansion is determined from the number of parameter lines just read in
		self._parameters['noslaters'] = counter - 1																				



	def getString(self):
		"""@return: string representation of the potential that can be read back in from a file"""
		lines = []
		lines.append(' &slater&')
		if (self._boundaries['cutoff'] != None):
			li = [str(self._parameters['cutoff']), str(self._boundaries['cutoff'][0]), str(self._boundaries['cutoff'][1])]
			lines.append('\t'.join(li))
		else:
			lines.append(str(self._parameters['cutoff']))
		for i in range(1, self._parameters['noslaters'] + 1):
			li = [str(self._parameters['n' + str(i)]), str(self._parameters['c' + str(i)])]

			if (self._boundaries['c' + str(i)] != None):
				li.append(str(self._boundaries['c' + str(i)][0]));li.append(str(self._boundaries['c' + str(i)][1]))
			
			li.append(str(self._parameters['z' + str(i)]))
			if (self._boundaries['z' + str(i)] != None):
				li.append(str(self._boundaries['z' + str(i)][0]));li.append(str(self._boundaries['z' + str(i)][1]))
			
			lines.append('\t'.join(li))
		
		return ' \n'.join(lines)


				
	def value(self, r):
		""" calculates the value of the slater expansion
		at a given distance r 
		@parameter r: float containing distance 
		@return: float containing function value """
		dummy = 0.0
		if (r <= self._parameters['cutoff']):
			for i in range(1, self._parameters['noslaters'] + 1):
				ci = self._parameters['c' + str(i)]
				zeta = self._parameters['z' + str(i)]
				normal = (2.0 * zeta) ** self._parameters['n' + str(i)] * numpy.sqrt(2.0 * zeta / self._factorial(int(2 * self._parameters['n' + str(i)])))
				dummy = dummy + ci * normal * r ** (self._parameters['n' + str(i)] - 1) * math.exp(-zeta * r)
		else:
			dummy = None
		return dummy
	

	# there is no built-in factorial in python 2.4
	def _factorial(self, x):
		""" @parameter x: integer fo calculate factorial from
		@return: float with factorial """
		res = 1
		for i in xrange(2, x + 1):
			res *= i
		return res


	def derivative(self, r):
		""" returns derivative of the slater expansion
		at a given distance r
		@parameter r: float containing distance 
		@return: float containing function derivative value """
	
		dummy = 0.0
		if (r <= self._parameters['cutoff']):
			# loop runs from 1 to no of slater functions
			for i in range(1, self._parameters['noslaters'] + 1):
				zeta = self._parameters['z' + str(i)]
				normal = (2.0 * zeta) ** self._parameters['n' + str(i)] * numpy.sqrt(2.0 * zeta / self._factorial(int(2.0 * self._parameters['n' + str(i)])))
				dummy = dummy + self._parameters['c' + str(i)] * normal * math.exp(-zeta * r) * ((self._parameters['n' + str(i)] - 1) * \
							r ** (self._parameters['n' + str(i)] - 2) - zeta * r ** (self._parameters['n' + str(i)] - 1))
		else:
			dummy = None

		return dummy
		
		
		
	def getCutoffs(self):
		""" @return: tuple with inner and outer cutoff """			
		return (0.0, self._parameters['cutoff'])
		


	def getMutables(self):
		""" @return: list with all parameters allowed to change """
		dummy = []
		if (self._boundaries['cutoff'] != None):
			dummy.append(self._parameters['cutoff'])

		for i in range(1, self._parameters['noslaters'] + 1):
			if (self._boundaries['c' + str(i)] != None):
				dummy.append(self._parameters['c' + str(i)])
			if (self._boundaries['z' + str(i)] != None):
				dummy.append(self._parameters['z' + str(i)])
		return dummy							



	def updateMutables(self, mutated):
		""" @parameter mutated: list with updated parameters to write back
		into self._parameters
		"""
		counter = 0
		# number of non None elements in boundaries
		for (dummy, value) in self._boundaries.items():
			if (value != None):
				counter = counter + 1
		if (len(mutated) != counter):
			raise CalcError("inconsistent number of parameters to update")
		counter = 0			
		if (self._boundaries['cutoff'] != None):
			self._parameters['cutoff'] = mutated[0]
			counter = counter + 1
			
		for i in range(1, self._parameters['noslaters'] + 1):
			if (self._boundaries['c' + str(i)] != None):
				self._parameters['c' + str(i)] = mutated[counter]
				counter = counter + 1
			if (self._boundaries['z' + str(i)] != None):
				self._parameters['z' + str(i)] = mutated[counter]
				counter = counter + 1
				

