##############################################################################
# velocityVerletOptimizer.py
# Part of comatsci computational materials science toolkit
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <janknaup@gmail.com>
#@organization: Bremen Center for Compuational Materials Science
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup  <janknaup@gmail.com>

from __future__ import print_function
from Optimizer import Optimizer,constants
import numpy

class velocityVerletOptimizer(Optimizer):
	"""minimize function by following gradients dynamically, tracking velocites.
	Has options to project velocities to onto new gradients and externally reset the 
	velocities"""
	
	
	
	def __init__(self, options):
		"""construt Velocity Verlet optimizer, minimize F(X) by iterating Velocity-Verlet algorithm.
		<em>All parameters inside options dictionary.</em>
		@param options dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF report convergence if abs of largest dF/dX component <=maxF,
		for Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		@param maxFRMS=-1 report convergence if RMS(dF/dX ) <= maxFRMS
		negative value means, never converge due to RMS derivative
		@param hardConverge=False only report convergence, if all convergence 
		criteria are met. This means, that convergence due to forces will never be met,
		if any criterion has a negative value.
		For Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		Default: False := report convergence, if any convergence criterion is met.
		@param stepSize (initial) step size for optimization
		@param adaptive=False if true, dynamically adapt stepsize depending 
		on velocity-force alignment. Alignment is checked before possible velocity
		projection, so that adaptive stepsize and projected velocities are compatible.
		@param stepAdaptFactor=1.618033988 factor to multiply or divide stepsize 
		by, when adapting. Must be >1.
		@param minStepSize=initial_stepSize/1000 do not reduce stepsize below this value
		@param maxStepSize=initial_stepSize*10 do not grow stepsive beyond this value
		@param growThreshold=0.9 if (dF' dot dF) > growThreshold, enlarge adaptive stepsize
		@param shrinkThreshold=0.5 if (dF' dot dF) < shrinkThreshold, shrink adaptive stepsize
		@param projectVelocities in each step, project the velocities into the new force
		direction, following the projected Velocity Verlet method
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print("Steepest Descent Optimizer: initializing.")
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: max derivative component: {0:f}.".format(self._maxF))
		
		self._maxFRMS=options.get("maxFRMS",-1.0)
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: max RMS derivative: {0:f}.".format(self._maxFRMS))
		
		self._hardConvergence=options.get("hardConvergence",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: using hard convergence: {0:s}.".format(str(self._maxF)))
		if self._hardConvergence and (self._maxF<=0 or self._maxFRMS<=0):
				if self._verbosity >= constants.VBL_SILENCE:
					print("Velocity Verlet Optimizer: Warning: hard convergence with negative force criterion. Will never converge due to force.")
		
		self._stepSize=options["stepSize"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: initial stepSize: {0:f}.".format(self._stepSize))
			
		self._adaptive=options.get("adaptive",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: using adaptive step size: {0:s}.".format(str(self._adaptive)))
			
		self._stepAdaptFactor=options.get("stepAdaptFactor",1.618033988)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print("Velocity Verlet Optimizer: step size adaption factor: {0:f}.".format(self._stepAdaptFactor))
	
		# save arithmetic operations: If stepsize is to remain > minStepSize, it must
		# be > minStepSize*stepAdaptFactor before scaling. Calculate the comparison value
		# one here
		self._minStepSize=options.get("minStepSize",self._stepSize/1000)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print("Velocity Verlet Optimizer: minimum adaptive step size: {0:f}.".format(self._minStepSize))
		self._minStepSize*=self._stepAdaptFactor
		
		# c.f. minStepSize
		self._maxStepSize=options.get("maxStepSize",self._stepSize*10)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print("Velocity Verlet Optimizer: maximum adaptive step size: {0:f}.".format(self._maxStepSize))
		self._minStepSize/=self._stepAdaptFactor

		self._growThreshold=options.get("growThreshold",0.9)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print("Velocity Verlet Optimizer: alignment threshold for stepsize growth: {0:f}.".format(self._growThreshold))
			
		self._shrinkThreshold=options.get("shrinkThreshold",0.5)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print("Velocity Verlet Optimizer: alignment threshold for stepsize shrink: {0:f}.".format(self._shrinkThreshold))
		
		self._projectVelocities=options.get("projectVelocities",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print("Velocity Verlet Optimizer: using projected Velocity Verlet: {0:s}.".format(str(self._projectVelocities)))
		
		#initialize velocities Vector and required shape of velocities array
		self._velocities=None
		#this will lead to acceptance of any array as inital velocities, once velocities have been
		#set for the first time, require the same shape for subsequent velocity arrays.
		self._velocityShape=None
		#initialize masses vector
		self._masses=None
		


	def setVelocities(self, dX):
		"""set internal velocities array, either for initial velocities or velocity 
		rescaling in dynamics simulations. Array shape is enforced to be flat 
		match inital array dimensions.
		<em>User must ensure that inital velocities are compatible to X vector dimensions</em>
		@param dX velocities array"""
		if self._velocities==None and dX.shape[2]==None:
			self._velocities=dX
			self._velocityShape=self._velocities.shape
			return
		elif dX.shape==self._velocities.shape:
			self._velocities=dX
			return
		else:
			raise "Velocity Verlet Optimizer: Velocity array shape mismatch."



	def getVelocities(self):
		"""return velocities array"""
		return self._velocities
	velocities=property(getVelocities,setVelocities)
	
	
	
	def setMasses(self, m):
		"""set internal masses array, only accept once, before first iteration step
		<em>User must ensure that inital masses are compatible to X vector dimensions</em>
		@param dX masses array"""
		if self._masses==None:
			self._masses=m
			return
		else:
			raise "Velocity Verlet Optimizer: Masses already set."



	def getMasses(self):
		"""return masses array"""
		return self._masses
	masses=property(getMasses,setMasses)
	
	
	
	def getStepSize(self):
		"""return step size used in last iteration"""
		return self._stepSize
	stepSize=property(getStepSize)
	
	
	
	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of Velocity Verlet
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X parameter vector to minimize <b>ignored</b>
		@param F function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF first derivative of function: dF/dX
		@param d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		RMS=numpy.sqrt(numpy.add.reduce(dF*dF)/dF.shape[0])
		maxF=max((-min(dF),max(dF)))
		if RMS < self._maxFRMS and maxF < self._maxF:
			self._convreason="Hard convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: convergence criterion met: {0:s}".format(self._convreason))
			return True
		elif not self._hardConvergence:
			if maxF < self._maxF:
				self._convreason="Max force component convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: convergence criterion met: {0:s}".format(self._convreason))
				return True
			if RMS < self._maxFRMS:
				self._convreason="Force RMS convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: convergence criterion met: {0:s}".format(self._convreason))
				return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print("Velocity Verlet Optimizer: convergence criterion met: {0:s}".format(self._convreason))
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print("Velocity Verlet Optimizer: not yet converged")
		return False

	
	
	def __adaptStepSize(self,dF):
		"""Adapt stepsize, using bounds, factors and thresholds defined upon 
		optimizer initialization"""
		#do nothing if no velocities are available or verlocity is zero
		if self._velocities==None or numpy.add.reduce(numpy.equal(self._velocities,0.0))==self._velocities.shape[0]:
			self._velocities=numpy.zeros(dF.shape,dtype=float)
			if self._verbosity >= constants.VBL_DEBUG2:
				print("Velocity Verlet Optimizer: No force-verlocity alignment info, skipping stepsize adaption")
			return
		#calculate derivative-velocity alignment
		align=numpy.dot(self._velocities,-dF)
		align/=numpy.sqrt(numpy.dot(self._velocities,self._velocities))
		align/=numpy.sqrt(numpy.dot(dF,dF))
		#grow, shrink or do nothing?
		if align <= self._shrinkThreshold:
			# only shrink, if stepsize is not below shrink minimum
			if self._stepSize > self._minStepSize:
				self._stepSize/=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: shrinking stepsize")
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: stepsize too small, skipping stepsize shrink")
		elif align >= self._growThreshold:
			if self._stepSize < self._maxStepSize:
				self._stepSize*=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: growing stepsize")
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print("Velocity Verlet Optimizer: stepsize too large, skipping stepsize grow")
		return



	def _velocitiesProject(self, dF):
		"""project velocities to derivative array to work as projected velocity verlet
		algorithm. velocities are stored in Class private array, hence no return value
		Set velocities to 0, if dX dot dF < 0.
		@param dF derivative array to project velocities onto
		"""
		#first calculate force direction
		forceDir=-dF/numpy.sqrt(numpy.dot(dF,dF))
		align=numpy.dot(forceDir, self._velocities)
		if align >=0:
			self._velocities=forceDir*align
		else:
			self._velocities=numpy.zeros(self._velocityShape,dtype=float)
	
	
	
	def _step(self,X,F,dF,d2F):
		"""perform one Velocity Verlet step, adapting stepsize and projecting velocities
		if required.
		c.f. base class documentation
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX <b>mandatory</b>
		@param d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#should we happen not to have forces or masses, set defaults now
		#this must happen here instead of in the initializer, since we do not know
		#the vector size at initialization time
		if self._velocities==None:
			self._velocities=numpy.zeros(X.shape,dtype=float)
			self._velocityShape=self._velocities.shape
		if self._masses==None:
			self._masses=numpy.ones(X.shape,dtype=float)
		#check if we have forces
		if dF==None:
			raise ValueError("derivative required for Velocity Verlet optimization")
		if self._adaptive:
			if self._verbosity >= constants.VBL_DEBUG2:
				print("Velocity Verlet Optimizer: checking stepsize")
			self.__adaptStepSize(dF)
		if self._projectVelocities:
			self._velocitiesProject(dF)
		#first update velocities
		self._velocities+=0.5*(-dF)*self._masses
		#then calculate new parameters
		newX=X+self._stepSize*self._velocities
		if self._verbosity >= constants.VBL_DEBUG2:
				print("Velocity Verlet Optimizer: step performed") 
		return newX



