## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# Optimizers.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function

#import verbosity levels from outside, if available, otherwise use own definitions
try:
	from .. import constants
except ImportError:
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
	@param element (fitness value, X) tuple"""
	return element[0]
	
	



class Optimizer:
	"""Function optimizer base class ("pure virtual") to define optimizer interface
	Optimizers minimize scalar function F(X), dependent on X, F(X) and optionally
	dF/dX and d2F/dX^2.
	Optimizers should not produce output below verosity level VBL_DEBUG1."""
	
	def __init__(self, options):
		"""optimizer base constructor
		always call in derived classes! <em>All parameters inside options dictionary.</em>
		@param options dictionary of optimizer options <b>contains all further parameters</b>
		@param verbosity=constants.VBL_QUIET output verbosity level, only warnings and fatal errors by default
		@param maxIterations Maximum number of interations until convergence due to iteration count is reported"""
		self._converged = False
		self._verbosity=options.get("verbosity",constants.VBL_QUIET)
		self._iterations=0
		self._maxIterations=options["maxIterations"]
		self._convreason="Not yet converged"
		if self._verbosity >= constants.VBL_DEBUG2:
			print("Optimizer basic initialization complete.")



	def _checkConvergence(self,X,F, dF, d2F):
		"""check convergence of optimization
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX, can be optional in implementation
		@param d2F second derivative of function d2F/dX2"""
		raise NotImplementedError
	
	
	
	def optStep(self, X, F, dF=None, d2F=None):
		"""Check convergence, if not converged, perform one single optimization step of F(X)
		Optimizer may provide multiple solutions. In that case, optStep returns one solution, all others are accessid via Optimizer.solutions property
		Implementations may make parameters dF,d2F optional.
		Do not reimplement, implement Optimizer._step() !
		Returns new vector X if not converged before optimization step, input vector X otherwise.
		@param X parameter vector to minimize
		@param F function value to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX, can be optional in implementation
		@param d2F second derivative of function d2F/dX2"""
		self._converged=self._checkConvergence(X,F,dF,d2F)
		if not self._converged:
			if self._verbosity >= constants.VBL_DEBUG1:
				print("Optimizer: not converged, performing iteration step")
			newX=self._step(X,F,dF,d2F)
			self._iterations+=1
		else:
			if self._verbosity >= constants.VBL_DEBUG1:
				print("Optimizer: converged, doing nothing")
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
		@param X parameter vector to minimize
		@param F function value to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX, can be optional in implementation
		@param d2F second derivative of function d2F/dX2"""
		raise NotImplementedError


