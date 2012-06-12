##############################################################################
# singleObjectiveMonteCarloOptimizer.py
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
import copy

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
		@param options dictionary of optimizer options <b>contains all further 
		parameters</b>, also see base class for global options!
		@param maxF report convergence if abs of largest dF/dX component <=maxF,
		for Dynamics simulation, set maxF and/or maxFRMS negative and set hardConvergence
		True
		@param mutator fuction with parameters mutator(options,X) returning 
		vector X suitable as argument of F(X)
		@param mutatorOptions dictionary of options to pass to the mutator function
		@param fitness fitness function F(X) with parameters F(options,X) returning 
		scalar fitness value
		@param fitnessOptions dictionary of options to pass to the fitness function
		"""
		if options["verbosity"] >= constants.VBL_DEBUG2:
			print("SOMC Optimizer: initializing.")
		#call base class constructor
		Optimizer.__init__(self,options)
		#digest opions, c.f. documentation
		self._maxF=options["maxF"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: max fitness: {0:f}.".format(self._maxF))

		self._mutator=options["mutator"]
##		if type(self._mutator)!="function":
##			raise "SOMC Optimizer: mutator function must be of type function"
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: mutator function stored")
		
		self._mutatorOptions=options["mutatorOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: mutator options stored")
		elif self._verbosity >= constants.VBL_DEBUG2:
			print("SOMC Optimizer: mutator options: {0:s}".format(str(self._mutatorOptions)))
		
		self._fitness=options["fitness"]
##		if type(self._fitness)!="function":
##			raise "SOMC Optimizer: fitness function must be of type function"
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: fitness function stored")
		
		self._fitnessOptions=options["fitnessOptions"]
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: fitness options stored")
		elif self._verbosity >= constants.VBL_DEBUG2:
			print("SOMC Optimizer: fitness options: {0:s}".format(str(self._fitnessOptions)))
		
		# initialize internal valiables
		# optimization specific
		self._lastFitness=None
		self._bestFitness=None
		self._bestX=None
		# statistics
		self._accepted=0
		self._mutations=0
		
		if self._verbosity >= constants.VBL_DEBUG1:
			print("SOMC Optimizer: initialized")
	
	
	
	def _checkConvergence(self,X,F,dF,d2F):
		"""check convergence of Monte-Carlo
		For internal use only, called in optStep method, convergence state is stored internally.
		External interface for convergence is the Optimizer.converged property
		Returns boolean, true if converged
		@param X parameter vector to minimize <b>ignored</b>
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX <b> ignored</b>
		@param d2F second derivative of function d2F/dX2 <b>ignored</b>"""
		#never report convergence in the first iteration, since no previous
		#fitness is stored when _checkConvergence is called...
		if self._iterations==0:
			return False
		#we have only one convergence criterion, which is maximum fitness
		if self._lastFitness < self._maxF:
			self._convreason="Max fitness convergence"
			if self._verbosity >= constants.VBL_DEBUG2:
				print("SOMC Optimizer: convergence criterion met: {0:s}".format(self._convreason))
			return True
		#also check number of iterations
		if self._iterations>=self._maxIterations:
			self._convreason="Maximum iterations reached"
			if self._verbosity >= constants.VBL_DEBUG2:
				print("SOMC Optimizer: convergence criterion met: {0:s}".format(self._convreason))
			return True
		if self._verbosity >= constants.VBL_DEBUG2:
			print("SOMC Optimizer: not yet converged")
		return False




	def _step(self,X,F,dF,d2F):
		"""perform one Monte-Carlo step, calculating fitness on the way
		c.f. base class documentation
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization<b>always ignored</b>
		@param dF first derivative of function: dF/dX <b>always ignored</b>
		@param d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
		#at first iteration, calculate fitness for input vector and return
		if self.iterations==0:
			if self._verbosity >= constants.VBL_DEBUG2:
				print("SOMC Optimizer: calculating fitness of initial X")
			self._bestFitness=self._fitness(self._fitnessOptions,X)
			self._lastFitness=copy.deepcopy(self._bestFitness)
			self._bestX=copy.deepcopy(X)
			self._lastSolutions=((self._bestFitness,self._bestX),)
			return X
		else:
		#otherwise mutate X, calculate new fitness and return the fitter X
			#mutate X
			if self._verbosity >= constants.VBL_DEBUG2:
				print("SOMC Optimizer: mutating X")
			newX=self._mutator(self._mutatorOptions,X)
			self._mutations+=1
			if self._verbosity >= constants.VBL_DEBUG2:
				print("SOMC Optimizer: calculating fitness of mutated X")
			self._lastFitness=self._fitness(self._fitnessOptions,newX)
			if self._lastFitness<self._bestFitness:
				self._bestFitness=copy.deepcopy(self._lastFitness)
				self._bestX=newX
				self._accepted+=1
				if self._verbosity >= constants.VBL_DEBUG2:
					print("SOMC Optimizer: mutated X is fitter, returning X'")
			else:
				if self._verbosity >= constants.VBL_DEBUG2:
					print("SOMC Optimizer: mutated X is not fitter, returning X")
				self._lastSolutions=copy.deepcopy(X)
			self._lastSolutions=((self._bestFitness,self._bestX),)
			return self._bestX


	
	def getAcceptanceRate(self):
		"""return total acceptance rate of trial soluations"""
		return float(self._accepted/float(self._mutations+0.000000000000000001))
	acceptanceRate=property(getAcceptanceRate)


