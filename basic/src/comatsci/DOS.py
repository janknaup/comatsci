## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# DOS.py - module to handle densities of states
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

import numpy.oldnumeric as num

import os,sys,math,re

import utils, constants

try:
	from scipy.signal import correlate
except ImportError:
	lfilter=None
	convolve=None
	
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
	
	spreadfunctions={
					"lorentz" :	lorentz,
					"gauss" :	gaussian
					}
	
	def  __init__(self):
		self.eigenValues=None
		self.fillings=None
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
		# initialize lists for temporary eigenvalue and fillings storage
		tempEigenValues=[[[]]]
		tempFillings=[[[]]]
		self.spins=1
		self.kpoints=1
		kindex=0
		spinindex=0
		# iterate through input file lines, ignoring line 1. Abort if line format indicates the beginning of a new k-point/spin section
		for i in range(1,len(eigenlines)):
			dummy=eigenlines[i].split()
			if len(dummy)>=4:
				try:
					if int(dummy[3])-1>spinindex:
						spinindex+=1
						kindex=0
						tempEigenValues.append([[]])
						tempFillings.append([[]])
					elif int(dummy[1])-1>kindex:
						kindex+=1
						tempEigenValues[spinindex].append([])
						tempFillings[spinindex].append([])
				except:
					print 'Error parsing line %d of file "%s" in band.out format. Abort.' % (i+1,filename)
					raise
				continue
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
				tempEigenValues[spinindex][kindex].append(float(dummy[0]))
			except:
				print 'Could not parse eigenvalue "%s" in line %d of file "%s".' %(dummy[0],i+1,filename)
				raise
			try:
				tempFillings[spinindex][kindex].append(float(dummy[1]))
			except:
				print 'Could not parse filling "%s" in line %d of file "%s".' %(dummy[1],i+1,filename)
				raise
			# end eigenvalue line parsing loop
		# store internal data
		self.eigenValues=num.array(tempEigenValues)
		self.fillings=num.array(tempFillings)
		self.spins=len(self.eigenValues)
		self.kpoints=len(self.eigenValues[0])
		# set dummy k-point weights
		self.kweights=num.ones((self.kpoints,), num.Float)/float(self.kpoints)
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
			print 'Could not open file "%s" in .EIG format. Abort.' % filename
			raise
		eigenlines=list(bandoutfile)
		bandoutfile.close()
		# we only read 1st k-point and first spin for the moment
		# initialize lists for temporary eigenvalue and fillings storage
		tempEigenValues=[[[]]]
		tempFillings=[[[]]]
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
				tempEigenValues[0][0].append(float(dummy[0]))
			except:
				print 'Could not parse eigenvalue "%s" in line %d of file "%s".' %(i+1,dummy[0],filename)
				raise
			# eig files contain no filling, so just store zero
			tempFillings[0][0].append(0.)
			# end eigenvalue line parsing loop
			# store internal data
			self.eigenValues=num.array(tempEigenValues,num.Float)
			self.fillings=num.array(tempFillings,num.Float)
			self.kweights=num.array([1.0,])
			# reset calculated properties, in case someone in reusing this instance
			self.reset_calculated()
			print self.eigenValues
			print self.fillings
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
		eigenShape=(eigenvalEntry.shape[2],eigenvalEntry.shape[1],eigenvalEntry.shape[0])
		# catch results.tag files, that do not contain eigenvalues:
		if eigenvalEntry==None:
			raise ValueError("No eigenvalues in tagged.out file '%s'." % filename)
		else:
			self.eigenValues=num.reshape(num.array(eigenvalEntry.value,num.Float),eigenShape)
			#self.eigenValues=num.reshape(num.array(eigenvalEntry.value,num.Float),eigenvalEntry.shape)
		self.eigenValues/=constants.EVOLT
		fillingsEntry=tagresults.getEntry("fillings")
		if fillingsEntry==None:
			self.fillings=num.ones(self.eigenValues.shape,num.Float)
		else:
			self.fillings=num.reshape(num.array(fillingsEntry.value,num.Float),eigenShape)
		if len(self.fillings) != len(self.eigenValues):
			raise ValueError("Number of fillings in tagged.out file '%s' does not match number of eigenvalues!" % filename)
		# set shape variables
		self.spins=len(self.eigenValues)
		self.kpoints=len(self.eigenValues[0])
		# set dummy k-point weights
		self.kweights=num.ones((self.kpoints,), num.Float)/float(self.kpoints)
		# reset calculated properties, in case someone in reusing this instance
		self.reset_calculated()
		# finished



	def readSiestaEIG(self, filename):
		"""Read eigenstates from SIESTA .EIG file
		@param filename: name of the file to read (can be .gz or .bz2 compressed)
		"""
		#open .EIG file, read all lines, concatenate and split into tokens
		try:
			eigFile=utils.compressedopen(filename,"r")
		except:
			print "could not open input file '%s'. abort." % (filename)
			raise
		eigLines=list(eigFile)
		eigFile.close()
		eigString=" ".join(eigLines)
		eigTokens=eigString.split()
		# parse header information
		self.fermiEnergy=float(eigTokens.pop(0))
		orbitalcount=int(eigTokens.pop(0))
		self.spins=int(eigTokens.pop(0))
		self.kpoints=int(eigTokens.pop(0))
		self.eigenValues=num.zeros((self.spins,self.kpoints,orbitalcount))
		self.fillings=num.zeros((self.spins,self.kpoints,orbitalcount))
		for k in range(self.kpoints):
			# drop k label
			eigTokens.pop(0)
			for s in range(self.spins):
				for i in range(orbitalcount):
					eigenValue=float(eigTokens.pop(0))
					self.eigenValues[s][k][i]=eigenValue
					if eigenValue<=self.fermiEnergy:
						self.fillings[s][k][i]=(2.0/float(self.spins))
					else:
						self.fillings[s][k][i]=0.0
		# set dummy k-point weights
		self.kweights=num.ones((self.kpoints,), num.Float)/float(self.kpoints)


	def hasEigenValues(self):
		"""check if Eigenvalues have been read
		@return: boolean, TRUE if eigenvalues have been read, FALSE otherwise
		"""
		return bool(self.eigenValues!=None)
	


	def readAimsOutput(self,filename):
		""" parse eigenspectrum from FHI aims standard output
		@type filename: string
		@param filename: name of the file containing the FHI aims console output  
		"""
		# prepare regurlar expressions
		pivotRE=re.compile("Writing Kohn-Sham eigenvalues.")
		# memory intensive implementation...
		infile=utils.compressedopen(filename, "r", autodetect=False)
		lines=list(infile)
		infile.close()
		# find _last_ occurence of Kohn-Sham Eigenstates marker, so iterate reverse from end
		# but use positive index counters for better readbility of forward iteration later
		i=len(lines)-1
		KSEpivot=None
		while(KSEpivot==None):
			if pivotRE.search(lines[i].strip())!=None:
				KSEpivot=i
			i-=1
		if KSEpivot==None: raise ValueError("No Kohn-SHam Eigenvalues found in aims output.")
		# Fermi level is output above eigenvalues list
		for i in range(KSEpivot-55,KSEpivot):
			if "(Fermi level)" in lines[i].strip():
				tokens=lines[i].strip().split()
				try:
					self.fermiEnergy=float(tokens[-1])
				except:
					print >> sys.stderr,"Failed to parse Fermi level in line %d of aims output file"%(i+1)
					raise
		# iterate below KSpivot to read eigenvalue lines.
		# Block of eigenvalues is enclosed in empty lines
		tempfillings=[]
		tempEigenvalues=[]
		i=KSEpivot+3
		while(lines[i].strip()!=""):
			tokens=lines[i].strip().split()
			try: 
				tempfillings.append(float(tokens[1]))
				tempEigenvalues.append(float(tokens[3]))
			except:
				print >> sys.stderr, "Failed to parse Kohn-Sham eigenvalues line in line %d of aims output file"%(i+1)
				raise
			i+=1
		# finished parsing, now store fillings and eigenvalues
		# FIXME: This implementation does not read k or spin resolved eigenstates
		self.kpoints=1
		self.spins=1
		self.kweights=num.array([1.0],num.Float)
		self.eigenValues=num.array([[tempEigenvalues]],num.Float)
		self.fillings=num.array([[tempfillings]],num.Float)
		
				



	###############################################################################
	# some useful properties
	###############################################################################
	
	def getMaxEigenValue(self):
		if self.eigenValues!=None:
			return num.amax(self.eigenValues)
		else:
			return None
	
	maxEigenValue=property(fget=getMaxEigenValue,fset=None,fdel=None,doc="""maximum Eigenvalue of eigenspectrum represented in this DOS instance""")
	
	def getMinEigenValue(self):
		if self.eigenValues!=None:
			return num.amin(self.eigenValues)
		else:
			return None
	
	minEigenValue=property(fget=getMinEigenValue,fset=None,fdel=None,doc="""minimum Eigenvalue of eigenspectrum represented in this DOS instance""")
	
	
	
	###############################################################################
	# Contiuous DOS representations
	###############################################################################
	
	def fastSpreadDOS(self, spreadfunction="lorentz",stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of spreading peaks
		NEW IMPLEMENTATION USING HISTOGRAM AND CORRELATE. MUCH FASTER BUT BINNING LEADS TO PEAK DISTORTIONS
		@type spreadfunction: string 
		@param spreadfunction: name of the mathematical function describing DOS peaks
			supported spread functions
			* lorentz - Cauchy-Lorentz distribution
			* gauss - Gaussian distribution
		@type stepwidth: float  
		@param stepwidth: step width of output array in eV (default 0.1)
		@type spread: float
		@param spread: peak spread of superposition, in eV (default 0.1)
		@type emin: float
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@type emax: float
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: tuple of Arrays, containing energies, and 2d arrays of total, occupied and unoccupied DOS per spin channel
		"""
		# check spreadfunction parameter
		if not spreadfunction in self.spreadfunctions:
			raise ValueError("Unknown spread function requested")
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
		emax+=((3*stepwidth)-(emax%stepwidth))  # add another stepwidtdh to emax to ensure highest eigenvalue is captured in binning
		# sanity check: step width must be smaller than peak width or weird things happen
		if (stepwidth>spread):
			raise ValueError("DOS sampling step width must be smaller than peak width")
		# sanity check: if emin >= emax, we are in trouble
		if emin >= emax:
			raise ValueError("emax not greater than emin")
		# initalize energies and DOS arrays
		energies=num.arange(emin,emax,stepwidth) # energies for DOS values
		print energies
		numsteps=len(energies)
		peakenergies=num.arange(-stepwidth*numsteps/2, stepwidth*numsteps/2, stepwidth) # energies for peak function
		# fill peak function array
		spreadpeak=num.fromfunction(lambda II: self.spreadfunctions[spreadfunction](x=peakenergies[II],s=spread,x0=0.0),(numsteps,))
		# initialize per psin channel DOS list
		spindoses=[]
		# define what filling is a full spin orbital (will probably break with non-colinear spin
		fullorbital=2/self.spins
		# iterate through spin channels
		for s in range(self.spins):
			# initialize per spin DOS arrays
			spreadDOS=num.zeros((numsteps,),num.Float)
			occDOS=num.zeros((numsteps,),num.Float)
			unoccDOS=num.zeros((numsteps,),num.Float)
			# iterate through k-points inside spin channel
			for k in range(self.kpoints):
				# bin eigenvalues by energy, use weighted histogram for occupation aware DOS
				doshist=num.histogram(self.eigenValues[s][k], energies)[0]
				occhist=num.histogram(self.eigenValues[s][k], energies, weights=self.fillings[s][k])[0]
				# convolve DOS histogram with peak function
				spreadDOStemp=(correlate(spreadpeak,doshist,"same")*fullorbital*self.kweights[k])[::-1]
				# accumulate spin-channel k-weighted DOS
				spreadDOS+=spreadDOStemp
				# convolve occupation aware DOS with peak function
				occDOStemp=(correlate(spreadpeak,occhist,"same")*self.kweights[k])[::-1]
				# accumulate
				occDOS+=occDOStemp
				# accumulate unoccupied DOS as difference total-occupied
				unoccDOS+=spreadDOStemp-occDOStemp #(correlate(spreadpeak,unocchist,"same")*self.kweights[k])[::-1]
			# store final array
			spindoses.append(num.array((spreadDOS,occDOS,unoccDOS),num.Float))
		# finished, return combined array
		return (num.array(energies),)+tuple(spindoses)
		
		
	
	def spreadDOS(self, spreadfunction="lorentz",stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of lorentz peaks
		@type spreadfunction: string 
		@param spreadfunction: name of the mathematical function describing DOS peaks
			supported spread functions
			* lorentz - Cauchy-Lorentz distribution
			* gauss - Gaussian distribution
		@type stepwidth: float  
		@param stepwidth: step width of output array in eV (default 0.1)
		@type spread: float
		@param spread: peak spread of superposition, in eV (default 0.1)
		@type emin: float
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@type emax: float
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: tuple of Arrays, containing energies, and 2d arrays of total, occupied and unoccupied DOS per spin channel
		"""
		# check spreadfunction parameter
		if not spreadfunction in self.spreadfunctions:
			raise ValueError("Unknown spread function requested")
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
		energies=num.arange(emin,emax,stepwidth)
		numsteps=len(energies)
		spindoses=[]
		fullorbital=2/self.spins
		for s in range(self.spins):
			spreadDOS=num.zeros((numsteps,),num.Float)
			occDOS=num.zeros((numsteps,),num.Float)
			unoccDOS=num.zeros((numsteps,),num.Float)
			# fill array by brute-force looping lorentz distributions for each eigenvalue
			for k in range(self.kpoints):
				for i in range(len(self.eigenValues[0][0])):
					for j in range(numsteps):
						rawdos=self.spreadfunctions[spreadfunction](energies[j],spread,self.eigenValues[s][k][i])*self.kweights[k]
						spreadDOS[j]+=rawdos
						occDOS[j]+=rawdos*self.fillings[s][k][i]
						unoccDOS[j]+=rawdos*(fullorbital-self.fillings[s][k][i])
			spindoses.append(num.array((spreadDOS,occDOS,unoccDOS),num.Float))
		# finished, return combined array
		return (num.array(energies),)+tuple(spindoses)
	
	
	
	
	def lorentzDOS(self,stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of lorentz peaks
		@param stepwidth: step width of output array in eV (default 0.1)
		@param spread: peak spread of superposition, in eV (default 0.1)
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: 2D-Array, containing energies in dimension 0 and DOS in dimension 1
		"""
		return self.spreadDOS("lorentz", stepwidth, spread, emin, emax)




	def gaussianDOS(self,stepwidth=0.1,spread=0.1,emin=None, emax=None):
		"""return de-discretized DOS using superposition of Gaussian peaks
		@param stepwidth: step width of output array in eV (default 0.1)
		@param spread: peak spread of superposition, in eV (default 0.1)
		@param emin: minimum energy of array to return, default: minimum eigenvalue (default None)
		@param emax: maximum energy of array to return, default: maximum eigenvalue (default None)
		@return: 2D-Array, containing energies in dimension 0 and DOS in dimension 1
		"""
		return self.spreadDOS("gauss", stepwidth, spread, emin, emax)




	def JDOS(self, spreadfunction="lorentz", stepwidth=0.1, spread=0.1):
		"""
		calculate joint density of states for given DOS object
		@type spreadfunction: string
		@param spreadfunction: string giving the spread function to use for the input DOS calculations
				Valid spread functions are:
				* "lorentz" - Cauchy-Lorentz Distribution
				* "gauss" - Gaussian distribution
		@type spread: float
		@param spread: spread width parameter for the spreading function
		@type stepwidth: float
		@param stepwidth: sampling step width for DOS and JDOS
		@rtype: 2D array of floats
		@return: 2D array containing energy shift in [0] and JDOS in [1]
		"""
		DOS=self.spreadDOS(spreadfunction, stepwidth, spread)
		# spin channel 0 is always present
		J=correlate(DOS[1][1],DOS[1][2],"full")
		# if DOS is spin polarized, add second spin channel JDOS to first
		if self.spins==2:
			J+=correlate(DOS[2][1],DOS[2][2],"full")
		nsteps=len(J)/2
		step=DOS[0][1]-DOS[0][0]
		stop=step*nsteps
		# construct output array
		JDOS=num.array((num.arange(0, stop, step),J[nsteps:0:-1]))
		# finisehd, return.
		return JDOS


class PDOS(DOS):
	"""Representation for a projected Density of states
	"""
	
	#define a set of expected AO lables for DFTB+ eigenvec.out files
	AOLabels=("S1","P1","P2","P3","D1","D2","D3","D4","D5","F1","F2","F3","F4","F5","F6","F7")
	
	
	def __init__(self):
		# generic simple multiple inheritance constructor
		for i in self.__class__.__bases__:
			i.__init__(self)
	
	
	
	def readEigenvecOut(self, filename="eigenvec.out"):
		"""reads basis function Mulliken populations per MO from DFTB+ eingencev.out file
		Eigenvalues must have already been loaded into the instance.
		@type filename: string
		@param filename: name of the input file, default "eigenvec.out"
		"""
		# check that a set of eigenvalues has been read before trying to read PDOS
		if not self.hasEigenValues():
			raise ValueError("Attempt to read eigenvectors without reading eigenvalues first.")
		# check if input file exists
		if not os.path.exists(filename):
			raise ValueError("Specified DFTB+ eigenvector file '%s' does not exist"%filename)
		# open compressed files transparently (but without automatic compressed file extension replacement)
		EVfile=utils.compressedopen(filename,autodetect=False)
		# read whole file into memory and close file object. May waste memory but mich easier to implement
		EVlines=EVfile.readlines()
		EVfile.close()
		# the first line contains junk (human readable header), 2nd line is empty
		#    check first line, if we find the expected junk
		if not EVlines[0][:60]=="Coefficients and Mulliken populations of the atomic orbitals":
			raise ValueError("did not found expected header in eigenvector file '%s'."%filename)
		# The number of eigenvectors should be the same as the number of eigenvalues
		# we have to determine the number of atoms and orbitals from the first eigenvalue block
		#   first check line 3 if it shows Eigenvector: 1 and Spin: 1
		print EVlines[2].split()
		if not len(EVlines[2].split()) >=3 or not(int(EVlines[2].split()[1])==1) or not (int(EVlines[2].split()[3][0])==1):
			raise ValueError("Did not find expected Eigenvector 1 Spin 1 component in line 3 of file '%s'"%filename)
		# initialize list of orbitals per atom, number of atoms
		orbitalsperAtom=[]
		tempcoefficients=[]
		tempmulliken=[]
		atomcount=0
		# start processing from fourth line
		currentline=3
		tempmulliken.append([])
		tempcoefficients.append([])
		orbitalindex=0
		while not (len(EVlines[currentline].split()) > 0 and EVlines[currentline].split()[0]=="Eigenvector:"):
			lineParts=EVlines[currentline].split()
			if len(lineParts)==0:
				if currentline!=3:
					orbitalsperAtom.append(orbitalindex)
				atomcount+=1
				orbitalindex=0
				currentline+=1
				continue
			else:
				if not lineParts[0].strip()==self.AOLabels[orbitalindex]:
					raise ValueError("Expected did not find expected orbital %s at line %d of file '%s'"%(self.AOLabels[orbitalindex],currentline+1,filename))
				tempcoefficients[-1].append(float(lineParts[1]))
				tempmulliken[-1].append(float(lineParts[2]))
				orbitalindex+=1
				currentline+=1
				continue
		# the above loop increased the atom count by 1 too much
		atomcount-=1
		orbitalsperatom=num.array(orbitalsperAtom,num.Float)
		# calculate the number of orbitals __before__ the index atom
		orbitalsum=orbitalsperatom.cumsum()-orbitalsperatom
		numOrbitals=int(orbitalsperatom.cumsum()[-1])
		# sanity-check if calculated number of orbitals is equal to number of eigenvalues
		if numOrbitals!=len(self.eigenValues):
			raise ValueError("Number of orbitals in eigenvector file does not equal number of eigenvalues from eigenvalue file.")
		# one eigenvector block contains one line per atomic orbital plus one blank line per atom plus two header lines
		blockLength=int(numOrbitals+atomcount+2)
		# we understand two possibilities for the eigenvector file length:
		#  either it is a one spin-component file or a two spin-component file which has double the number of eigenvalue blocks
		if len(EVlines)==2+blockLength*numOrbitals:		# the file has 2 header lines before the eigenvector blocks
			numSpins=1
		elif len(EVlines)==2+2*blockLength*numOrbitals:
			numSpins=2
		else:
			raise ValueError("Unexpected number of lines in file '%s' file is probably malformed."%filename)
		# allocate data arrays
		self.orbitalCoeficients=num.zeros(numSpins*numOrbitals**2,num.Float)
		self.orbitalMulliken=num.zeros(numSpins*numOrbitals**2,num.Float)
		# store the number of orbitals and orbitals per atom as member variables
		self.numOrbitals=numOrbitals
		self.orbitalsPerAtom=orbitalsperatom
		del tempcoefficients
		del tempmulliken
		# read all eigenvector blocks (re-reading the first block, this is easier than catching the case of two spin components)
		# iterate eigenvalues
		for EVindex in range(self.numOrbitals):
			# iterate spin components per eigenvalue
			for spinIndex in range(numSpins):
				blockBase=2+blockLength*(EVindex+spinIndex)
				# check for expected eigenvalue block header (file indices count from 1)
				if not len(EVlines[blockBase].split()) >=3 or not(int(EVlines[blockBase].split()[1])==EVindex+1) or not (int(EVlines[blockBase].split()[3][0])==spinIndex+1):
					raise ValueError("Did not find expected Eigenvector %d Spin %d component in file '%s'"%(EVindex+1,spinIndex+1,filename))
				# iterate atoms per spin eigenvalue
				for atomIndex in range(atomcount):
					# iterate orbitals per atom
					for aoIndex in range(int(self.orbitalsPerAtom[atomIndex])):
						# calculate coefficient index
						coeffIndex=(EVindex+spinIndex)*numOrbitals+orbitalsum[atomIndex]+aoIndex
						#split the current line
						lineIndex=blockBase+int(orbitalsum[atomIndex])+2+atomIndex+aoIndex
						lineParts=EVlines[lineIndex].split()
						if not lineParts[0].strip()==self.AOLabels[aoIndex]:
							raise ValueError("Did not find expected atomic orbital '%s' in line %d of egenvector file '%s'"%(self.AOLabels[aoIndex],lineIndex+1,filename))
						else:
							try:
								self.orbitalCoeficients[coeffIndex]=float(lineParts[1])
								self.orbitalMulliken[coeffIndex]=float(lineParts[2])
							except:
								raise ValueError("failed to parse line %d of eigenvector file '%s'"%(lineIndex+1,filename))
								raise
		# done





	