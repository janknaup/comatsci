## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# DOS.py - module to handle densities of states
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

import numpy.oldnumeric as num

import gzip,bz2,os,shutil,sys,math

import utils, constants


###############################################################################
# some mathematical helper functions
###############################################################################

def lorentz(x,s=1.,x0=0.):
	"""Cauchy-Lorentz distribution function
	@param x: independent variable
	@param s: half-width parameter (default 1.0)
	@param x0: peak position (default 0.0)
	"""
	return (1./constants.PI)*(s/((s**2)+((x-x0)**2)))



def gaussian(x,s=1.,x0=0.):
	"""Gaussian (or normal) distribution function
	@param x: independent variable
	@param s: half-width parameter (default 1.0)
	@param x0: peak position (default 0.0)
	"""
	return (1./(s*math.sqrt(2*constants.PI)))*math.exp(-((x-x0)**2)/(2*s**2))



###############################################################################
# The class to handle eigenspektrum-derived DOSes
###############################################################################

class DOS:
	"""Class representing a density of states"""
	
	
	def  __init__(self):
		pass



	def reset_calculated(self):
		"""delete all calculated properties"""
		pass




	def readBandOut(self, filename, fakeFillings=False):
		"""Read eigenstates from DFTB+ band.out file
		@param filename: name of the file to read (can be .gz or .bz2 compressed)
		@param fakeFillings: if true, silently set missing filling values to 0.0 (default False)
		"""
		# try to open the input file and read all lines to memory
		try:
			bandoutfile=utils.compressedopen(filename,"r")
		except:
			print 'Could not open file "%s" in band.out format. Abort.' % filename
			raise
		eigenlines=list(bandoutfile)
		bandoutfile.close()
		# we only read 1st k-point and first spin for the moment
		# TODO: extend to multiple k-points and spins
		# initialize lists for temporary eigenvalue and fillings storage
		tempEigenValues=[]
		tempFillings=[]
		# iterate through input file lines, ignoring line 1. Abort if line format indicates the beginning of a new k-point/spin section
		for i in range(1,len(eigenlines)):
			dummy=eigenlines[i].split()
			if len(dummy)==4:
				# new section starts, abort reading.
				break
			elif len(dummy)==0:
			  # ignore empty lines
			  continue
			elif len(dummy)==1 and fakeFillings:
				# dirty trick: on caller's request, fake filling values of zero, if filling column is missing
				dummy.append("0.0")
			elif len(dummy)!=2:
				# we expect one energy and one filling, if we find anything more or less, abort.
				print 'Unexpected entry in line %d of "%s", cannot parse line "%s". Abort.' %(i+1,filename,eigenlines[i])
				sys.exit(1)
			# try to convert energy and filling to floats
			try:
				tempEigenValues.append(float(dummy[0]))
			except:
				print 'Could not parse eigenvalue "%s" in line %d of file "%s".' %(dummy[0],i+1,filename)
				raise
			try:
				tempFillings.append(float(dummy[1]))
			except:
				print 'Could not parse filling "%s" in line %d of file "%s".' %(dummy[1],i+1,filename)
				raise
			# end eigenvalue line parsing loop
			# store internal data
			self.eigenValues=num.array(tempEigenValues,num.Float)
			self.fillings=num.array(tempFillings,num.Float)
			# reset calculated properties, in case someone in reusing this instance
			self.reset_calculated()
			# finished


	def readDFTBEigFile(self, filename):
		"""Read eigenstates from (old) DFTB .eig file
		@param filename: name of the file to read (can be .gz or .bz2 compressed)
		"""
		# try to open the input file and read all lines to memory
		try:
			bandoutfile=utils.compressedopen(filename,"r")
		except:
			print 'Could not open file "%s" in band.out format. Abort.' % filename
			raise
		eigenlines=list(bandoutfile)
		bandoutfile.close()
		# we only read 1st k-point and first spin for the moment
		# initialize lists for temporary eigenvalue and fillings storage
		tempEigenValues=[]
		tempFillings=[]
		# iterate through input file lines, ignoring line 1. Abort if line format indicates the beginning of a new k-point/spin section
		for i in range(1,len(eigenlines)):
			# ignore comment lines
			if i.strip()[0] in ("#",):
				continue
			dummy=eigenlines[i].split()
			if len(dummy)==0:
			  # ignore empty lines
			  continue
			elif len(dummy)!=1:
				# we expect exactly one energy, if we find anything more , abort.
				print 'Unexpected entry in line %d of "%s", cannot parse line "%s". Abort.' %(i+1,filename,eigenlines[i])
				sys.exit(1)
			# try to convert energy and filling to floats
			try:
				tempEigenValues.append(float(dummy[0]))
			except:
				print 'Could not parse eigenvalue "%s" in line %d of file "%s".' %(i+1,dummy[0],filename)
				raise
			# eig files contain to filling, so just store zero
			tempFillings.append(0.)
			# end eigenvalue line parsing loop
			# store internal data
			self.eigenValues=num.array(tempEigenValues,num.Float)
			self.fillings=num.array(tempFillings,num.Float)
			# reset calculated properties, in case someone in reusing this instance
			self.reset_calculated()
			# finished


	def readTaggedOut(self, filename):
		"""Read eigenstates from DFTB+ tagged.out file
		@param filename: name of the file to read (can be .gz or .bz2 compressed)
		"""
		# try to open the input file and read all lines to memory
		try:
			resultsfile=utils.compressedopen(filename,"r")
		except:
			print 'Could not open file "%s" in band.out format. Abort.' % filename
			raise
		# now parse the tagged.out file using Balint's code
		tagresults=utils.TaggedCollection(utils.ResultParser(resultsfile).entries)
		# we don't need the file objecy anymore
		resultsfile.close()
		# extract eigenvalues and fillings from the results collection and store them
		eigenvalEntry=tagresults.getEntry("eigenvalues")
		# catch results.tag files, that do not contain eigenvalues:
		if eigenvalEntry==None:
			raise ValueError("No eigenvalues in tagged.out file '%s'." % filename)
		else:
			self.eigenValues=num.array(eigenvalEntry.value[0],num.Float)
		fillingsEntry=tagresults.getEntry("fillings")
		if eigenvalEntry==None:
			self.fillings=num.ones(fillingsEntry.shape[0],num.Float)
		else:
			self.fillings=num.array(fillingsEntry.value[0],num.Float)
		if len(self.fillings) != len(self.eigenValues):
			raise ValueError("Number of fillings in tagged.out file '%s' does not match number of eigenvalues!" % filename)
		# reset calculated properties, in case someone in reusing this instance
		self.reset_calculated()
		# finished



	def readSiestaEIG(self, filename):
		"""Read eigenstates from SIESTA .EIG file
		@param filename: name of the file to read (can be .gz or .bz2 compressed)
		"""
		pass


	###############################################################################
	# some useful properties
	###############################################################################
	
	def getMaxEigenValue(self):
		if self.eigenValues!=None:
			return max(self.eigenValues)
		else:
			return None
	
	maxEigenValue=property(fget=getMaxEigenValue,fset=None,fdel=None,doc="""maximum Eigenvalue of eigenspectrum represented in this DOS instance""")
	
	def getMinEigenValue(self):
		if self.eigenValues!=None:
			return min(self.eigenValues)
		else:
			return None
	
	minEigenValue=property(fget=getMinEigenValue,fset=None,fdel=None,doc="""minimum Eigenvalue of eigenspectrum represented in this DOS instance""")
	
	
	
	###############################################################################
	# Contiuous DOS representations
	###############################################################################
	
	def lorentzDOS(self,stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of lorentz peaks
		@param stepwidth: step width of output array in eV (default 0.1)
		@param spread: peak spread of superposition, in eV (default 0.1)
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: 2D-Array, containing energies in dimension 0 and DOS in dimension 1
		"""
		# sanity check: if we do not have eigenvalues, calculating a DOS cannot work:
		if self.eigenValues==None:
			raise RuntimeError("No eigenvalues to calculate DOS from")
		# if minimum energy is not specified, use minimum eigenvalue
		if emin==None:
			emin=self.minEigenValue
		# if maximum energy is not specified, use maximum eigenvalue
		if emax==None:
			emax=self.maxEigenValue
		# shift emin down and emax up to be multiples of stepwidth
		# this is necessary to ensure compatibility of DOSes of different eigenspectra
		emin-=emin % stepwidth
		emax+=stepwidth-emax%stepwidth
		# sanity check: if emin >= emax, we are in trouble
		if emin >= emax:
			raise ValueError("emax not greater than emin")
		# initalize energies DOS arrays
		energies=num.arrayrange(emin,emax,stepwidth)
		numsteps=len(energies)
		lorentzDOS=num.zeros((numsteps,),num.Float)
		# fill array by brute-force looping lorentz distributions for each eigenvalue
		for i in self.eigenValues:
			for j in range(numsteps):
			  lorentzDOS[j]+=lorentz(energies[j],spread,i)
		# finished, return combined array
		return num.array((energies,lorentzDOS),num.Float)



	def gaussianDOS(self,stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of Gaussian peaks
		@param stepwidth: step width of output array in eV (default 0.1)
		@param spread: peak spread of superposition, in eV (default 0.1)
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: 2D-Array, containing energies in dimension 0 and DOS in dimension 1
		"""
		# sanity check: if we do not have eigenvalues, calculating a DOS cannot work:
		if self.eigenValues==None:
			raise RuntimeError("No eigenvalues to calculate DOS from")
		# if minimum energy is not specified, use minimum eigenvalue
		if emin==None:
			emin=self.minEigenValue
		# if maximum energy is not specified, use maximum eigenvalue
		if emax==None:
			emax=self.maxEigenValue
		# shift emin down and emax up to be multiples of stepwidth
		# this is necessary to ensure compatibility of DOSes of different eigenspectra
		emin-=emin % stepwidth
		emax+=stepwidth-emax%stepwidth
		# sanity check: if emin >= emax, we are in trouble
		if emin >= emax:
			raise ValueError("emax not greater than emin")
		# initalize energies DOS arrays
		energies=num.arrayrange(emin,emax,stepwidth)
		numsteps=len(energies)
		gaussDOS=num.zeros((numsteps,),num.Float)
		# fill array by brute-force looping lorentz distributions for each eigenvalue
		for i in self.eigenValues:
			for j in range(numsteps):
			  gaussDOS[j]+=gaussian(energies[j],spread,i)
		# finished, return combined array
		return num.array((energies,gaussDOS),num.Float)

