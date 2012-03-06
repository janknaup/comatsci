##############################################################################
# muellerBrownCalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from comatsci.Calculators.Calculator import Calculator,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED#,CALCSTATUS_ERROR,CALCSTATUS_DISABLED

#from comatsci.Calculators.CalcError import CalcError
import comatsci.constants as constants
#import sys

import numpy


class muellerBrownCalc(Calculator):
	"""calculate energies and forces from a Mueller-Brown potential"""



	def __init__(self, verbosity=0):
		"""Initialize Mueller Brown potential calculator
		@param verbosity: c.f. base class (default 1)
		@type laplace: boolean
		@param laplace: calculate laplacian of the potential (delfault False)
"""
		Calculator.__init__(self,verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print "initializing Mueller-Brown calculator"
		self._status=CALCSTATUS_READY


	# We don't need to implement most of the template functions
	# since we don't call an external caclulator program
	# We must reimplement start(), finreready(), errorreready()
	# and runfg() though

	def start(self,Geometry, steplabel, charge=0):
		"""Calculte energies and forces.
		@param Geometry: c.f. base class
		@param steplabel: c.f. base class
		@param charge: ignored (default 0)
		prior to calculation"""
		# set potential parameters as in Mueller, Brown: Theor. Chim. Acta. 53 pp. 75ff.
		A=[-200.,-100.,-170.,15]
		a=[-1.,-1.,-6.5,0.7]
		b=[0.,0.,11.,0.6]
		c=[-10.,-10.,-6.5,0.7]
		x0=[1.,0.,-0.5,-1.]
		y0=[0.,0.5,1.5,1.]
		# set status flag to "running"
		self._status=CALCSTATUS_RUNNING
		# initialize forces array and laplace and energy variables
		tempforces=numpy.zeros((Geometry.Atomcount,3),dtype=float)
		tempenergy=0.0
		templaplace=0.0
		# iterate through atoms
		for i in range(Geometry.Atomcount):
			# iterate through Mueller-Brown potential terms
			for j in range(4):
				#save the exponential value, as we need it for the forces
				expoterm=A[j]*numpy.exp(a[j]*(Geometry.Geometry[i][0]-x0[j])**2 + b[j]*(Geometry.Geometry[i][0]-x0[j])*(Geometry.Geometry[i][1]-y0[j])+c[j]*(Geometry.Geometry[i][1]-y0[j])**2)
				# store energy contribution
				tempenergy+=expoterm
				# calculate and store froce contributions, store derivative prefactor for use in laplacian
				vx=(2*a[j]*(Geometry.Geometry[i][0]-x0[j])+b[j]*(Geometry.Geometry[i][1]-y0[j]))
				vy=(2*c[j]*(Geometry.Geometry[i][1]-y0[j])+b[j]*(Geometry.Geometry[i][0]-x0[j]))
				tempforces[i][0]-=expoterm*vx
				tempforces[i][1]-=expoterm*vy
				# calculate laplacian contribution
				templaplace+=expoterm*(vx*2*a[j]+vy*2*c[j])
		# finish iterations
		self.etot=tempenergy
		self.gradients=tempforces
		self.laplace=templaplace
		self._status=CALCSTATUS_FINISHED



	def runfg(self,Geometry, steplabel, charge=0, Erep=None):
		"""directly calculates energy and forces, see start()"""
		self.start(Geometry, steplabel, charge)



	def getplace(self):
		"""@return: laplacian vector from last calculation, if laplacian calculation was enabled upon initialization
		"""
		if not self._calcLaplace:
			raise RuntimeError("Attempt to retrieve laplacian from Mueller-Brown calculator that was not initialized to calculate it.")
		else:
			return self.laplace
	

