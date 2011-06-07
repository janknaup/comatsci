##############################################################################
# newtonRaphsonOptimizer.py
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


from Optimizer import *


class newtonRaphsonOptimizer(Optimizer):
	"""minimize function using Newton-Raphon algorithm
	"""
	
	
	
	def __init__(self, options):
		"""construt Newton-Raphson optimizer, minimize F(X) by iterating:
		X'=X-F(X)*(dF/|dF|**2) where F(X) must be scalar function of vector X
		and dF must be vector of partial derivatives dF/dX_i.
		<em>All parameters inside options dictionary.</em>
		@param options dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF report convergence if abs of largest dF/dX component <=maxF,
		@param maxFRMS=-1 report convergence if RMS(dF/dX ) <= maxFRMS
		negative value means, never converge due to RMS derivative
		@param hardConverge=False only report convergence, if all convergence 
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
		@param X parameter vector to minimize <b>ignored</b>
		@param F function to minimize, F can be vector for multi-objective optimization <b>ignored</b>
		@param dF first derivative of function: dF/dX
		@param d2F second derivative of function d2F/dX2 <b>ignored</b>"""
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
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX <b>mandatory</b>
		@param d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		delta=dF/num.dot(dF,dF)
		delta*=F
		return X-delta
		
		
