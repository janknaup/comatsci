##############################################################################
# SteepestDescentOptimizer.py
# Part of comatsci computational materials science toolkit
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <janknaup@gmail.com>
#@organization: EPFL - Ecole polytechnique federale de Lausanne
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup  <janknaup@gmail.com>


from Optimizer import Optimizer,constants


class steepestDescentOptimizer(Optimizer):
	"""Optimize X so that F(X) becomes zero by following dF/dX"""
	
	def __init__(self, options):
		"""construt steepest descent optimizer, minimize F(X) by iterating:
		X'=X-stepSize*(dF/dX)
		<em>All parameters inside options dictionary.</em>
		@param options dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF report convergence if abs of largest dF/dX component <=maxF
		@param maxFRMS=-1 report convergence if RMS(dF/dX ) <= maxFRMS
		negative value means, never converge due to RMS derivative
		@param hardConverge=False only report convergence, if all convergence 
		criteria are met. Must never be combined with negative force convergence criteria.
		Default: report convergence, if any convergence criterion is met.
		@param stepSize (initial) step size for optimization
		@param adaptive=False if true, dynamically adapt stepsize depending 
		on curvature
		@param constantDisplacement=False if true, normalize dF/dX to unity 
		in each iteration. Strongly recommended in conjunction with constantDisplacement
		@param stepAdaptFactor=1.618033988 factor to multiply or divide stepsize 
		by, when adapting. Must be >1.
		@param minStepSize=initial_stepSize/1000 do not reduce stepsize below this value
		@param maxStepSize=initial_stepSize*10 do not grow stepsive beyond this value
		@param growThreshold=0.9 if (dF' dot dF) > growThreshold, enlarge adaptive stepsize
		@param shrinkThreshold=0.5 if (dF' dot dF) < shrinkThreshold, shrink adaptive stepsize
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
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX <b>mandatory</b>
		@param d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
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
		@param X parameter vector to minimize <b>ignored</b>
		@param F function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF first derivative of function: dF/dX
		@param d2F second derivative of function d2F/dX2 <b>ignored</b>"""
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
	
