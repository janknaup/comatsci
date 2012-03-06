##############################################################################
# polynom.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <janknaup@gmail.com>
#@organization: Bremen Center for Compuational Materials Science
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup	<janknaup@gmail.com>

from potentialFunction import potentialFunction
from comatsci.Calculators.CalcError import *


class polynomFunction(potentialFunction):
	""" class for polynom represenation of a pairpotential """

	def __init__(self,parameters):
		""" polynom constructor
		@parameter parameters: dictionary of function parameters """

		potentialFunction.__init__(self,parameters)



	def parseString(self,inString):
		""" read polynom coefficients from inString string
		update self._parameters
		@parameter inString: string containing function parameters, separated by \n """
		inlines=inString.split("\n")
		""" first line should contain &polynom& """
		dummy=inlines[0].split()[0].strip()
		if not (dummy=='&polynom&'):
			raise CalcError(" format of polynom coefficient file not recognized ")

		dummy=''
		counter=0 # counts the number of lines read-in as holding parameters
		# remember: the first line should contain the potential function identifier field, e.g. &polynom&

		for i in inlines[1:len(inlines)]:
		# ignore common comment markers and blank lines 
			if not (i==''):
				firstchar=i.lstrip()[0]
				if not (firstchar=="#" or firstchar==";" or firstchar=="!" or firstchar=='*'):
					dummy=i.split()
					counter=counter+1
					# the first parameters containing line holds the cutoff
					if (counter==1):
						self._parameters['cutoff']=float(dummy[0])
						if (len(dummy)==3):
									# if there are two more numbers in the line
									# the cutoff is assumed to be subject to optimization
									# within the given boundaries
							self._boundaries['cutoff']=(float(dummy[1]),float(dummy[2]))
						elif (len(dummy)==1):
									# no boundaries given: cutoff will be held fixed
									# during optimization
							self._boundaries['cutoff']=None
						else:
							raise CalcError(" incorrect no of variable parameter boundaries ")
					else:
						self._parameters['c'+str(counter-2)]=float(dummy[0])
						if (len(dummy)==3):
							self._boundaries['c'+str(counter-2)]=(float(dummy[1]),float(dummy[2]))
						elif (len(dummy)==1):
							self._boundaries['c'+str(counter-2)]=None
						else:
							raise CalcError(" incorrect no of variable parameter boundaries ")
		# maximum power in the polynom is determined from the number of parameters just read in
		self._parameters['order']=counter-2																				


	def getString(self):
		"""@return string representation of the potential that can be read back in from a file"""
		lines=[]
		lines.append(' &polynom&')
		if (self._boundaries['cutoff'] != None):
			li=[str(self._parameters['cutoff']),str(self._boundaries['cutoff'][0]),str(self._boundaries['cutoff'][1])]
			lines.append('\t'.join(li))
		else:
			lines.append(str(self._parameters['cutoff']))
		for i in range(self._parameters['order']+1):
			if (self._boundaries['c'+str(i)] != None):
				li=[str(self._parameters['c'+str(i)]),str(self._boundaries['c'+str(i)][0]),str(self._boundaries['c'+str(i)][1])]
				lines.append('\t'.join(li))
			else:
				lines.append(str(self._parameters['c'+str(i)]))
		
		return ' \n'.join(lines)
				
	def value(self,r):
		""" calculates the value of the polynom
		at a given distance r 
		@parameter r: float containing distance 
		@return float containing function value """
		dummy=0.0
		if (r <= self._parameters['cutoff']):
			for i in range(self._parameters['order']+1):
				dummy=dummy+self._parameters['c'+str(i)]*(self._parameters['cutoff']-r)**i
		return dummy


	def derivative(self,r):
		""" returns derivative of the polynom
		at a given distance r
		@parameter r: float containing distance 
		@return float containing function derivative value """
	
		dummy=0.0
		if (r <= self._parameters['cutoff']):
			for i in range(self._parameters['order']+1):
				dummy=dummy-i*self._parameters['c'+str(i)]*(self._parameters['cutoff']-r)**(i-1)
		return dummy
		
		
		
	def getCutoffs(self):
		""" @return tuple with inner and outer cutoff """			
		return (1E-10,self._parameters['cutoff']-1E-10)
		


	def getMutables(self):
		""" @return list with all parameters allowed to change """
		dummy=[]
		if (self._boundaries['cutoff'] != None):
			dummy.append(self._parameters['cutoff'])

		for i in range(self._parameters['order']+1):
			if (self._boundaries['c'+str(i)] !=None):
				dummy.append(self._parameters['c'+str(i)])
		return dummy							



	def updateMutables(self,mutated):
		""" @parameter mutated: list with updated parameters to write back
		into self._parameters
		"""
		counter=0
		# number of non None elements in boundaries
		for (dummy,value) in self._boundaries.items():
			if (value != None):
				counter=counter+1
		if (len(mutated) != counter):
			raise CalcError("inconsistent number of parameters to update")
		counter=0			
		if (self._boundaries['cutoff'] != None):
			self._parameters['cutoff']=mutated[0]
			counter=counter+1
			
		for i in range(self._parameters['order']+1):
			if (self._boundaries['c'+str(i)] !=None):
				self._parameters['c'+str(i)]=mutated[counter]
				counter=counter+1
