##############################################################################
# singleObjectiveGeneticOptimizer.py
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
		@param combiner function that combines two X vectors, returning the child,
		must have interface combiner(options, X1, X2)
		@param combinerOptions dictionary of options to pass to the combiner
		function
		@param populationSize=20 number of solution candidates in population
		@param breederCount=5 number of breeders to select from population
		@param keepFitterElders=False if true, keep all members of the parent generation
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
		@param X parameter vector to minimize <b>ignored</b>
		@param F function to minimize, F can be vector for multi-objective optimization
		@param dF first derivative of function: dF/dX <b> ignored</b>
		@param d2F second derivative of function d2F/dX2 <b>ignored</b>"""
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
		@param X parameter vector to minimize
		@param F function to minimize, F can be vector for multi-objective optimization<b>always ignored</b>
		@param dF first derivative of function: dF/dX <b>always ignored</b>
		@param d2F second derivative of function d2F/dX2 <b>always ignored</b>"""
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
		@param X initial input vector to generate population from"""
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


