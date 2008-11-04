## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# Optimizers.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

import numpy.oldnumeric as num
import math, copy, random

#import verbosity levels from outside, if available, otherwise use own definitions
try:
	import constants
except:
	class constants:
		"""default Verbosity level definitions"""
		# Verbosity levels:
		# only fatal errors
		VBL_SILENCE=-1
		# fatal errors and warnings
		VBL_QUIET=0
		# VBL_QUIET+status and progress reports
		VBL_NORMAL=1
		#VBL_NORMAL+extended status reports
		VBL_TALKY=2
		#VBL_TALKY+basic debug information
		VBL_DEBUG1=3
		#VBL_DEBUG1+extended debug information
		VBL_DEBUG2=4



def populationKey(element):
	"""helper function to sort (fitness,X) tuples by fitness
	@param element: (fitness value, X) tuple"""
	return element[0]
	
	



class Optimizer:
	"""Function optimizer base class ("pure virtual") to define optimizer interface
	Optimizers minimize scalar function F(X), dependent on X, F(X) and optionally
	dF/dX and d2F/dX^2.
	Optimizers should not produce output below verosity level VBL_DEBUG1."""
	
	def __init__(self, options):
		"""optimizer base constructor
		always call in derived classes! <em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further parameters</b>
		@param: verbositty=constants.VBL_QUIET output verbosity level, only warnings and fatal errors by default
		@param maxIterations: Maximum number of interations until convergence due to iteration count is reported"""
		self._converged = False
		self._verbosity=options.get("verbosity",constants.VBL_QUIET)
		self._iterations=0
		self._maxIterations=options["maxIterations"]
		self._convreason="Not yet converged"
		if self._verbosity >= constants.VBL_DEBUG2:
			print "Optimizer basic initialization complete."



	def _checkConvergence(self,X,F, dF, d2F):
		"""check convergence of optimization
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX, can be optional in implementation
		@param: d2F second derivative of function d2F/dX2"""
		pass
	
	
	
	def optStep(self, X, F, dF=None, d2F=None):
		"""Check convergence, if not converged, perform one single optimization step of F(X)
		Optimizer may provide multiple solutions. In that case, optStep returns one solution, all others are accessid via Optimizer.solutions property
		Implementations may make parameters dF,d2F optional.
		Do not reimplement, implement Optimizer._step() !
		Returns new vector X if not converged before optimization step, input vector X otherwise.
		@param X: parameter vector to minimize
		@param F: function value to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX, can be optional in implementation
		@param: d2F second derivative of function d2F/dX2"""
		self._converged=self._checkConvergence(X,F,dF,d2F)
		if not self._converged:
			if self._verbosity >= constants.VBL_DEBUG1:
				print "Optimizer: not converged, performing iteration step"
			newX=self._step(X,F,dF,d2F)
			self._iterations+=1
		else:
			if self._verbosity >= constants.VBL_DEBUG1:
				print "Optimizer: converged, doing nothing"
			newX=X
		return newX
	
	
	
	def getConverged(self):
		"""return true if convergence criteria defined in options are met"""
		return self._converged
	converged=property(getConverged)
	
	
	
	def getSolutions(self):
		"""return last optimized parameter array"""
		return self._lastSolutions
	solutions=property(getSolutions)
	
	
	
	def getIterations(self):
		"""return number of iterations performed so far"""
		return self._iterations
	iterations=property(getIterations)
	
	
	
	def getMaxIterations(self):
		"""return maximum number of iterations"""
		return self._maxIterations
	maxIterations=property(getMaxIterations)
	
	
	
	def getConvReason(self):
		"""return convergence reason string"""
		return self._convreason
	convReason=property(getConvReason)
	
	
	
	def _step(self,X,F,dF,d2F):
		"""really perform one single optimization step of F(X)
		Optimizer may provide multiple solutions. In that case, optStep returns one solution, all others are accessid via Optimizer.solutions property
		Implementations may make parameters dF,d2F optional.
		<b>always reimplement this!</b>
		Returns new vector X.
		@param X: parameter vector to minimize
		@param F: function value to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX, can be optional in implementation
		@param: d2F second derivative of function d2F/dX2"""
		pass




class steepestDescentOptimizer(Optimizer):
	"""Optimize X so that F(X) becomes zero by following dF/dX"""
	
	def __init__(self, options):
		"""construt steepest descent optimizer, minimize F(X) by iterating:
		X'=X-stepSize*(dF/dX)
		<em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF: report convergence if abs of largest dF/dX component <=maxF
		@param: maxFRMS=-1 report convergence if RMS(dF/dX ) < maxFRMS (default )
		negative value means, never converge due to RMS derivative
		@param hardConverge: only report convergence, if all convergence  (default False)
		criteria are met. Must never be combined with negative force convergence criteria.
		Default: report convergence, if any convergence criterion is met.
		@param stepSize: (initial) step size for optimization
		@param adaptive: if true, dynamically adapt stepsize depending  (default False)
		on curvature
		@param constantDisplacement: if true, normalize dF/dX to unity  (default False)
		in each iteration. Strongly recommended in conjunction with constantDisplacement
		@param stepAdaptFactor: factor to multiply or divide stepsize  (default 1.618033988)
		by, when adapting. Must be >1.
		@param: minStepSize=initial_stepSize/1000 do not reduce stepsize below this value
		@param: maxStepSize=initial_stepSize*10 do not grow stepsive beyond this value
		@param growThreshold: if (dF' dot dF) > growThreshold, enlarge adaptive stepsize (default 0.9)
		@param shrinkThreshold: if (dF' dot dF) < shrinkThreshold, shrink adaptive stepsize (default 0.5)
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print "Steepest Descent Optimizer: initializing."
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: max Force component: %f." % (self._maxF)
			
		self._maxFRMS=options.get("maxFRMS",-1.0)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: max RMS Force: %f." % (self._maxFRMS)
			
		self._hardConvergence=options.get("hardConvergence",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: using hard convergence: %s." % (str(self._hardConvergence))
		if self._hardConvergence and (maxF<=0 or maxFRMS<=0):
##			raise "Hard convergence demands all positive derivative convergence criteria"
				if self._verbosity >= constants.VBL_SILENCE:
					print "Steepest Descent Optimizer: Warning: hard convergence with negative froce criterion. Will never converge due to force."
			
		self._stepSize=options["stepSize"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: initial step size: %f." % (self._stepSize)
			
		self._adaptive=options.get("adaptive",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: adaptive stepsize: %s." % (str(self._adaptive))
			
		self._constantDisplacement=options.get("constantDisplacement",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Steepest Descent Optimizer: derivative indepentent displacement: %s." % (str(self._constantDisplacement))
			
		self._stepAdaptFactor=options.get("stepAdaptFactor",1.618033988)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Steepest Descent Optimizer: stepsize adaption factor: %f." % (self._stepAdaptFactor)
		if self._adaptive and self._stepAdaptFactor<=1:
			raise "Steepest descent optimizer: stepsiza adaption factor must be > 1"
			
		# save arithmetic operations: If stepsize is to remain > minStepSize, it must
		# be > minStepSize*stepAdaptFactor before scaling. Calculate the comparison value
		# one here
		self._minStepSize=options.get("minStepSize",self._stepSize/1000)*self._stepAdaptFactor
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Steepest Descent Optimizer: minimum adaptive stepsize: %f." % (self._minStepSize/self._stepAdaptFactor)
		
		# c.f. minStepSize above
		self._maxStepSize=options.get("maxStepSize",self._stepSize*10)/self._stepAdaptFactor
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Steepest Descent Optimizer: maximuim adaptive stepsize: %f." % (self._maxStepSize*self._stepAdaptFactor)
			
		self._growThreshold=options.get("growThreshold",0.9)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Steepest Descent Optimizer: alignment threshold for stepsize growth: %f." % (self._growThreshold)
			
		self._shrinkThreshold=options.get("shrinkThreshold",0.5)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Steepest Descent Optimizer: alignment threshold for stepsize shrink: %f." % (self._shrinkThreshold)
			
		#initialize history needed for stepsize adaption, convergence checking etc.
		self._oldForces=None
		if self._verbosity >= constants.VBL_DEBUG2:
			print "Steepest Descent Optimizer: Initialized."



	def __adaptStepSize(self,dF):
		"""Adapt stepsize, using bounds, factors and thresholds defined upon 
		optimizer initialization"""
		#do nothing if no old forces are stored
		if self._oldForces==None:
			self._oldForces=copy.deepcopy(dF)
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Steepest Descent Optimizer: No force alignment info, skipping stepsize adaption"
			return
		#calculate derivative alignment
		align=num.dot(self._oldForces,dF)
		align/=math.sqrt(num.dot(self._oldForces,self._oldForces))
		align/=math.sqrt(num.dot(dF,dF))
		#grow, shrink or do nothing?
		if align <= self._shrinkThreshold:
			# only shrink, if stepsize is not below shrink minimum
			if self._stepSize > self._minStepSize:
				self._stepSize/=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: shrinking stepsize"
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: stepsize too small, skipping stepsize shrink"
		elif align >= self._growThreshold:
			if self._stepSize < self._maxStepSize:
				self._stepSize*=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: growing stepsize"
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: stepsize too large, skipping stepsize grow"
		self._oldForces=copy.deepcopy(dF)
		return
	
	
	
	def _step(self,X,F,dF,d2F):
		"""perform one steepest descent step, adapting stepsize if required
		c.f. base class documentation
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX <b>mandatory</b>
		@param: d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#check if we have forces
		if dF==None:
			raise "derivative required for steepest descent optimization"
		if self._constantDisplacement:
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Steepest Descent Optimizer: normalizing derivative vector"
			dF/=math.sqrt(num.dot(dF,dF))
		if self._adaptive:
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Steepest Descent Optimizer: checking stepsize"
			self.__adaptStepSize(dF)
		newX=X-dF*self._stepSize
		if self._verbosity >= constants.VBL_DEBUG2:
				print "Steepest Descent Optimizer: steepest descent step performed" 
		return newX



	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of steepest descent
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X: parameter vector to minimize <b>ignored</b>
		@param F: function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF: first derivative of function: dF/dX
		@param: d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		RMS=math.sqrt(num.add.reduce(dF*dF)/dF.shape[0])
		maxF=max((-min(dF),max(dF)))
		if RMS < self._maxFRMS and maxF < self._maxF:
			self._convreason="Hard convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		elif not self._hardConvergence:
			if maxF < self._maxF:
				self._convreason="Max force component convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
			if RMS < self._maxFRMS:
				self._convreason="Force RMS convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Steepest Descent Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Steepest Descent Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print "Steepest Descent Optimizer: not yet converged"
		return False
	
	
	def getStepSize(self):
		"""return step size used in last iteration"""
		return self._stepSize
	stepSize=property(getStepSize)
	
	




class velocityVerletOptimizer(Optimizer):
	"""minimize function by following gradients dynamically, tracking velocites.
	Has options to project velocities to onto new gradients and externally reset the 
	velocities"""
	
	
	
	def __init__(self, options):
		"""construt Velocity Verlet optimizer, minimize F(X) by iterating Velocity-Verlet algorithm.
		<em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF: report convergence if abs of largest dF/dX component <=maxF,
		for Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		@param: maxFRMS=-1 report convergence if RMS(dF/dX ) < maxFRMS (default )
		negative value means, never converge due to RMS derivative
		@param hardConverge: only report convergence, if all convergence  (default False)
		criteria are met. This means, that convergence due to forces will never be met,
		if any criterion has a negative value.
		For Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		Default: False := report convergence, if any convergence criterion is met.
		@param stepSize: (initial) step size for optimization
		@param adaptive: if true, dynamically adapt stepsize depending  (default False)
		on velocity-force alignment. Alignment is checked before possible velocity
		projection, so that adaptive stepsize and projected velocities are compatible.
		@param stepAdaptFactor: factor to multiply or divide stepsize  (default 1.618033988)
		by, when adapting. Must be >1.
		@param: minStepSize=initial_stepSize/1000 do not reduce stepsize below this value
		@param: maxStepSize=initial_stepSize*10 do not grow stepsive beyond this value
		@param growThreshold: if (dF' dot dF) > growThreshold, enlarge adaptive stepsize (default 0.9)
		@param shrinkThreshold: if (dF' dot dF) < shrinkThreshold, shrink adaptive stepsize (default 0.5)
		@param projectVelocities: in each step, project the velocities into the new force
		direction, following the projected Velocity Verlet method
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print "Steepest Descent Optimizer: initializing."
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: max derivative component: %f." % (self._maxF)
		
		self._maxFRMS=options.get("maxFRMS",-1.0)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: max RMS derivative: %f." % (self._maxFRMS)
		
		self._hardConvergence=options.get("hardConvergence",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: using hard convergence: %s." % (str(self._maxF))
		if self._hardConvergence and (maxF<=0 or maxFRMS<=0):
				if self._verbosity >= constants.VBL_SILENCE:
					print "Velocity Verlet Optimizer: Warning: hard convergence with negative froce criterion. Will never converge due to force."
		
		self._stepSize=options["stepSize"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: initial stepSize: %f." % (self._stepSize)
			
		self._adaptive=options.get("adaptive",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: using adaptive step size: %s." % (str(self._adaptive))
			
		self._stepAdaptFactor=options.get("stepAdaptFactor",1.618033988)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Velocity Verlet Optimizer: step size adaption factor: %f." % (self._stepAdaptFactor)
	
		# save arithmetic operations: If stepsize is to remain > minStepSize, it must
		# be > minStepSize*stepAdaptFactor before scaling. Calculate the comparison value
		# one here
		self._minStepSize=options.get("minStepSize",self._stepSize/1000)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Velocity Verlet Optimizer: minimum adaptive step size: %f." % (self._minStepSize)
		self._minStepSize*=self._stepAdaptFactor
		
		# c.f. minStepSize
		self._maxStepSize=options.get("maxStepSize",self._stepSize*10)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Velocity Verlet Optimizer: maximum adaptive step size: %f." % (self._maxStepSize)
		self._minStepSize/=self._stepAdaptFactor

		self._growThreshold=options.get("growThreshold",0.9)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Velocity Verlet Optimizer: alignment threshold for stepsize growth: %f." % (self._growThreshold)
			
		self._shrinkThreshold=options.get("shrinkThreshold",0.5)
		if self._verbosity >= constants.VBL_DEBUG1 and self._adaptive:
			print "Velocity Verlet Optimizer: alignment threshold for stepsize shrink: %f." % (self._shrinkThreshold)
		
		self._projectVelocities=options.get("projectVelocities",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Velocity Verlet Optimizer: using projected Velocity Verlet: %s." % (str(self._projectVelocities))
		
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
		@param dX: velocities array"""
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
		@param dX: masses array"""
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
		@param X: parameter vector to minimize <b>ignored</b>
		@param F: function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF: first derivative of function: dF/dX
		@param: d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		RMS=math.sqrt(num.add.reduce(dF*dF)/dF.shape[0])
		maxF=max((-min(dF),max(dF)))
		if RMS < self._maxFRMS and maxF < self._maxF:
			self._convreason="Hard convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		elif not self._hardConvergence:
			if maxF < self._maxF:
				self._convreason="Max force component convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
			if RMS < self._maxFRMS:
				self._convreason="Force RMS convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Velocity Verlet Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print "Velocity Verlet Optimizer: not yet converged"
		return False

	
	
	def __adaptStepSize(self,dF):
		"""Adapt stepsize, using bounds, factors and thresholds defined upon 
		optimizer initialization"""
		#do nothing if no velocities are available or verlocity is zero
		if self._velocities==None or num.add.reduce(num.equal(self._velocities,0.0))==self._velocities.shape[0]:
			self._velocities=num.zeros(dF.shape,num.Float)
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Velocity Verlet Optimizer: No force-verlocity alignment info, skipping stepsize adaption"
			return
		#calculate derivative-velocity alignment
		align=num.dot(self._velocities,-dF)
		align/=math.sqrt(num.dot(self._velocities,self._velocities))
		align/=math.sqrt(num.dot(dF,dF))
		#grow, shrink or do nothing?
		if align <= self._shrinkThreshold:
			# only shrink, if stepsize is not below shrink minimum
			if self._stepSize > self._minStepSize:
				self._stepSize/=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: shrinking stepsize"
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: stepsize too small, skipping stepsize shrink"
		elif align >= self._growThreshold:
			if self._stepSize < self._maxStepSize:
				self._stepSize*=self._stepAdaptFactor
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: growing stepsize"
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Velocity Verlet Optimizer: stepsize too large, skipping stepsize grow"
		return



	def _velocitiesProject(self, dF):
		"""project velocities to derivative array to work as projected velocity verlet
		algorithm. velocities are stored in Class private array, hence no return value
		Set velocities to 0, if dX dot dF < 0.
		@param dF: derivative array to project velocities onto
		"""
		#first calculate force direction
		forceDir=-dF/math.sqrt(num.dot(dF,dF))
		align=num.dot(forceDir, self._velocities)
		if align >=0:
			self._velocities=forceDir*align
		else:
			self._velocities=num.zeros(self._velocityShape,num.Float)
	
	
	
	def _step(self,X,F,dF,d2F):
		"""perform one Velocity Verlet step, adapting stepsize and projecting velocities
		if required.
		c.f. base class documentation
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX <b>mandatory</b>
		@param: d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#should we happen not to have forces or masses, set defaults now
		#this must happen here instead of in the initializer, since we do not know
		#the vector size at initialization time
		if self._velocities==None:
			self._velocities=num.zeros(X.shape,num.Float)
			self._velocityShape=self._velocities.shape
		if self._masses==None:
			self._masses=num.ones(X.shape,num.Float)
		#check if we have forces
		if dF==None:
			raise "derivative required for Velocity Verlet optimization"
		if self._adaptive:
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Velocity Verlet Optimizer: checking stepsize"
			self.__adaptStepSize(dF)
		if self._projectVelocities:
			self._velocitiesProject(dF)
		#first update velocities
		self._velocities+=0.5*(-dF)*self._masses
		#then calculate new parameters
		newX=X+self._stepSize*self._velocities
		if self._verbosity >= constants.VBL_DEBUG2:
				print "Velocity Verlet Optimizer: step performed" 
		return newX






class singleObjectiveMonteCarloOptimizer(Optimizer):
	"""minimize function using random walk algorithm
	c.f. base class documentation"""
	
	
	
	def __init__(self, options):
		"""SOMC optimizer, minimize F(X) by randomly modifying X. F(X) is to be understood 
		as a single objective fitness function of X, where X is fitter than X' if F(X)<F(X').<br>
		F(X) must be passed to the optimizer as a callable object, with the actual interface:
		F(options,X), where options is a dictionary of options to the fitness function.<br>
		<em>All parameters F, dF, d2F passed to public methods of SOMC are dummies</em><br>
		<em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF: report convergence if abs of largest dF/dX component <=maxF,
		for Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		@param mutator: fuction with parameters mutator(options,X) returning 
		vector X suitable as argument of F(X)
		@param mutatorOptions: dictionary of options to pass to the mutator function
		@param fitness: fitness function F(X) with parameters F(options,X) returning 
		scalar fitness value
		@param fitnessOptions: dictionary of options to pass to the fitness function
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print "SOMC Optimizer: initializing."
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: max fitness: %f." % (self._maxF)

		self._mutator=options["mutator"]
##		if type(self._mutator)!="function":
##			raise "SOMC Optimizer: mutator function must be of type function"
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: mutator function stored"
		
		self._mutatorOptions=options["mutatorOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: mutator options stored"
		elif self._verbosity >= constants.VBL_DEBUG2:
			print "SOMC Optimizer: mutator options: %s" % (str(self._mutatorOptions))
		
		self._fitness=options["fitness"]
##		if type(self._fitness)!="function":
##			raise "SOMC Optimizer: fitness function must be of type function"
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: fitness function stored"
		
		self._fitnessOptions=options["fitnessOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: fitness options stored"
		elif self._verbosity >= constants.VBL_DEBUG2:
			print "SOMC Optimizer: fitness options: %s" % (str(self._fitnessOptions))
		
		# initialize internal valiables
		# optimization specific
		self._lastFitness=None
		self._bestFitness=None
		self._bestX=None
		# statistics
		self._accepted=0
		
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOMC Optimizer: initialized"
	
	
	
	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of Monte-Carlo
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X: parameter vector to minimize <b>ignored</b>
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX <b> ignored</b>
		@param: d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		#never report convergence in the first iteration, since no previous
		#fitness is stored when _checkConvergence is called...
		if self._iterations==0:
			return False
		#we have only one convergence criterion, which is maximum fitness
		if self._lastFitness < self._maxF:
			self._convreason="Max fitness convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOMC Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOMC Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print "SOMC Optimizer: not yet converged"
		return False




	def _step(self,X,F,dF,d2F):
		"""perform one Monte-Carlo step, calculating fitness on the way
		c.f. base class documentation
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization<b>always ignored</b>
		@param dF: first derivative of function: dF/dX <b>always ignored</b>
		@param: d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#at first iteration, calculate fitness for input vector and return
		if self.iterations==0:
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOMC Optimizer: calculating fitness of initial X"
			self._bestFitness=self._fitness(self._fitnessOptions,X)
			self._lastFitness=copy.deepcopy(self._bestFitness)
			self._bestX=copy.deepcopy(X)
			return X
		else:
		#otherwise mutate X, calculate new fitness and return the fitter X
			#mutate X
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOMC Optimizer: mutating X"
			newX=self._mutator(self._mutatorOptions,X)
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOMC Optimizer: calculating fitness of mutated X"
			self._lastFitness=self._fitness(self._fitnessOptions,newX)
			if self._lastFitness<self._bestFitness:
				self._bestFitness=copy.deepcopy(self._lastFitness)
				if self._verbosity >= constants.VBL_DEBUG2:
					print "SOMC Optimizer: mutated X is fitter, returning X'"
				return newX
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print "SOMC Optimizer: mutated X is not fitter, returning X"
				return X


	
	def getAcceptanceRate(self):
		"""return total acceptance rate of trial soluations"""
		return float(self._accepted/float(self.iterations))
	acceptanceRate=property(getAcceptanceRate)





class singleObjectiveGeneticOptimizer(Optimizer):
	"""minimize function by breeding and selecting within a population of solutuion candidates
	"""
	
	def __init__(self, options):
		"""SOGA optimizer: minimize F(X) by randomly modifying and interbreeding X within a 
		population of solution candidates. F(X) is to be understood 
		as a single objective fitness function of X, where X is fitter than X' if F(X)<F(X').<br>
		F(X) must be passed to the optimizer as a callable object, with the actual interface:
		F(options,X), where options is a dictionary of options to the fitness function.<br>
		<em>All parameters F, dF, d2F passed to public methods of SOMC are dummies</em><br>
		<em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF: report convergence if abs of largest dF/dX component <=maxF,
		for Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		@param mutator: fuction with parameters mutator(options,X) returning 
		vector X suitable as argument of F(X)
		@param mutatorOptions: dictionary of options to pass to the mutator function
		@param fitness: fitness function F(X) with parameters F(options,X) returning 
		scalar fitness value
		@param fitnessOptions: dictionary of options to pass to the fitness function
		@param combiner: function that combines two X vectors, returning the child,
		must have interface combiner(options, X1, X2)
		@param combinerOptions: dictionary of options to pass to the combiner
		function
		@param populationSize: number of solution candidates in population (default 20)
		@param breederCount: number of breeders to select from population (default 5)
		@param keepFitterElders: if true, keep all members of the parent generation (default False)
		that are fitter than the children, otherwise keep exactly the breeders, replace
		all other specimens
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print "SOGA Optimizer: initializing."
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: max fitness: %f." % (self._maxF)

		self._mutator=options["mutator"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: mutator function stored"
		
		self._mutatorOptions=options["mutatorOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: mutator options stored"
		elif self._verbosity >= constants.VBL_DEBUG2:
			print "SOGA Optimizer: mutator options: %s" % (str(self._mutatorOptions))
		
		self._fitness=options["fitness"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: fitness function stored"
		
		self._fitnessOptions=options["fitnessOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: fitness options stored"
		elif self._verbosity >= constants.VBL_DEBUG2:
			print "SOGA Optimizer: fitness options: %s" % (str(self._fitnessOptions))
		
		self._combiner=options["combiner"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: combiner stored"
		
		self._combinerOptions=options["combinerOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: combiner options stored"
		elif self._verbosity >= constants.VBL_DEBUG2:
			print "SOGA Optimizer: combiner options: %s" % (str(self._combinerOptions))
		
		self._populationSize=options.get("populationSize",20)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: population size: %d." % (self._populationSize)
		
		self._breederCount=options.get("breederCount",5)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: breeder Count: %d." % (self._breederCount)
		
		self._keepFitterElders=options.get("keepFitterElders",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: keeping fitter elders: %s." % (str(self._keepFitterElders))
		
		# initialize internal valiables
		self._lastPopulation=None
		self._bestPopulation=None
		
		# statistics
		self._mutations=0
		self._accepted=0
		
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: initialized"
			
	

	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of SOGA
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X: parameter vector to minimize <b>ignored</b>
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX <b> ignored</b>
		@param: d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		#never report convergence in the first iteration, since no previous
		#fitness is stored when _checkConvergence is called...
		if self._iterations==0:
			return False
		#we have only one convergence criterion, which is maximum fitness
		if self._lastPopulation[0][0] < self._maxF:
			self._convreason="Max fitness convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOGA Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOGA Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print "SOGA Optimizer: not yet converged"
		return False




	def _step(self,X,F,dF,d2F):
		"""perform one Monte-Carlo step, calculating fitness on the way
		c.f. base class documentation
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization<b>always ignored</b>
		@param dF: first derivative of function: dF/dX <b>always ignored</b>
		@param: d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#at first iteration, generate initial population and return
		if self.iterations==0:
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOGA Optimizer: generating inital population"
			self.__initialPopulate(X)
			return self._bestPopulation[0][1]
		else:
		#otherwise combine, mutate, and evaluate required number of offspring
			if self._verbosity >= constants.VBL_DEBUG2:
				print "SOGA Optimizer: breeding new generation"
			newPopulation=copy.deepcopy(self._lastPopulation[0:self._breederCount])
			for i in range (self._breederCount,self._populationSize):
				#combine
				newX=self._combiner(self._combinerOptions,
						random.choice(self._lastPopulation[0:self._breederCount])[1],
						random.choice(self._lastPopulation[0:self._breederCount])[1])
				#mutate, calc fitness and store
				newX=self._mutator(self._mutatorOptions,newX)
				self._mutations+=1
				newPopulation.append(tuple((self._fitness(self._fitnessOptions,newX),newX)))
			#now sort the new population by fitness
			newPopulation.sort(key=populationKey)
			#if we want to kill off all non-breeder elders, we are finished,
			#otherwise we have to merge the old non-breeders into the new
			#population, re-sort and slice to the desired population size
			if not self._keepFitterElders:
				self._lastPopulation=newPopulation
			else:
				newPopulation+=self._lastPopulation[self._breederCount:]
				newPopulation.sort(key=populationKey)
				self._lastPopulation=newPopulation[:self._populationSize]
			return self._lastPopulation[0][1]


	
	def getAcceptanceRate(self):
		"""return total acceptance rate of trial soluations"""
		return float(self._accepted/float(self.iterations))
	acceptanceRate=property(getAcceptanceRate)


	
	def __initialPopulate(self, X):
		"""generate initial population for genetic optimization
		@param X: initial input vector to generate population from"""
		if self._verbosity >= constants.VBL_DEBUG1:
			print "SOGA Optimizer: generating initial population and fitnesses"
		#We store the population as a list of (F(X),X) tuples
		#calculate the fitness of the initial X here and store in population list
		self._lastPopulation=[(self._fitness(self._fitnessOptions,X),X),]
		# now generate populationSize-1 mutants andcalculate their fitnesses
		for i in range(1, self._populationSize):
			newX=self._mutator(self._mutatorOptions,X)
			self._mutations+=1
			self._lastPopulation.append(tuple((self._fitness(self._fitnessOptions,newX),newX)))
		# now sort population by fitness
		self._lastPopulation.sort(key=populationKey)
		#at last, copy inital population into best population array
		self._bestPopulation=copy.deepcopy(self._lastPopulation)



class newtonRaphsonOptimizer(Optimizer):
	"""minimize function using Newton-Raphon algorithm
	"""
	
	
	
	def __init__(self, options):
		"""construt Newton-Raphson optimizer, minimize F(X) by iterating:
		X'=X-F(X)*(dF/|dF|**2) where F(X) must be scalar function of vector X
		and dF must be vector of partial derivatives dF/dX_i.
		<em>All parameters inside options dictionary.</em>
		@param options: dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF: report convergence if abs of largest dF/dX component <=maxF,
		@param: maxFRMS=-1 report convergence if RMS(dF/dX ) < maxFRMS (default )
		negative value means, never converge due to RMS derivative
		@param hardConverge: only report convergence, if all convergence  (default False)
		criteria are met. This means, that convergence due to forces will never be met,
		if any criterion has a negative value.
		Default: False := report convergence, if any convergence criterion is met.
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print "Newton-Raphson Optimizer: initializing."
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Newton-Raphson Optimizer: max derivative component: %f." % (self._maxF)
		
		self._maxFRMS=options.get("maxFRMS",-1.0)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Newton-Raphson Optimizer: max RMS derivative: %f." % (self._maxFRMS)
		
		self._hardConvergence=options.get("hardConvergence",False)
		if self._verbosity >= constants.VBL_DEBUG1:
			print "Newton-Raphson Optimizer: using hard convergence: %s." % (str(self._maxF))
		if self._hardConvergence and (maxF<=0 or maxFRMS<=0):
			raise "Newton-Raphson Optimizer: Warning: hard convergence with negative froce criterion. Will never converge due to force."





	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of Newton-Raphson optimization
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X: parameter vector to minimize <b>ignored</b>
		@param F: function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF: first derivative of function: dF/dX
		@param: d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		RMS=math.sqrt(num.add.reduce(dF*dF)/dF.shape[0])
		maxF=max((-min(dF),max(dF)))
		if RMS < self._maxFRMS and maxF < self._maxF:
			self._convreason="Hard convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
					print "Newton-Raphson Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		elif not self._hardConvergence:
			if maxF < self._maxF:
				self._convreason="Max force component convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Newton-Raphson Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
			if RMS < self._maxFRMS:
				self._convreason="Force RMS convergence"
				if self._verbosity >= constants.VBL_DEBUG2:
					print "Newton-Raphson Optimizer: convergence criterion met: %s" %(self._convreason)
				return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print "Newton-Raphson Optimizer: convergence criterion met: %s" %(self._convreason)
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print "Newton-Raphson Optimizer: not yet converged"
		return False



	def _step(self, X, F, dF, d2F):
		"""perform one Newton-Raphson step
		c.f. base class documentation
		@param X: parameter vector to minimize
		@param F: function to minimize, F can be vector for multi-objective optimization
		@param dF: first derivative of function: dF/dX <b>mandatory</b>
		@param: d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		delta=dF/num.dot(dF,dF)
		delta*=F
		return X-delta
		
		
