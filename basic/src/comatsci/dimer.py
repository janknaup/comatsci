## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# Dimer.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
import numpy

#rename geometry to cure naming Awkwardnesses
from .geometry.geometry import Geometry,FMG_DTD,GeometryError

import constants
import utils

import math,copy,os,sys
import xml.dom.minidom

DIMER_DTD=FMG_DTD+"""
<!ELEMENT DeltaR (#PCDATA)>
<!ATTLIST DeltaR
	lunit (ang|au) "ang"
>
<!ELEMENT NoGradInRot EMPTY>
<!ELEMENT curvature (#PCDATA)>
<!ELEMENT E0 (#PCDATA)>
<!ATTLIST E0
	eunit (eV|au|H) "au"
>
<!ELEMENT E1 (#PCDATA)>
<!ATTLIST E1
	eunit (eV|au|H) "au"
>
<!ELEMENT E2 (#PCDATA)>
<!ATTLIST E2
	eunit (eV|au|H) "au"
>
<!ELEMENT f0 (#PCDATA)>
<!ELEMENT f1 (#PCDATA)>
<!ELEMENT f2 (#PCDATA)>
<!ELEMENT fN (#PCDATA)>
<!ELEMENT Dimer (Geometry,DeltaR,NoGradInRot?,(E0,E1)?,E2?,(f0,f1)?,f2?,fN?,curvature?)>
"""

class Dimer(Geometry):
	"""Represents a dimer of two atomic configurations for transistion state searches
	using the advanced dimer method by Heyden et. al.: J. Phys. Chem. 123, 224101 (2005)"""



	def __init__(self, options=None):
		"""Inititalize Geometry object
		@keyword options: dictionary containing the options, use all defaults if None. <b>parameters below are passed in options dictionary!</b><br /> (default None)
		@keyword Mode:	Geometry mode, known modes are:
		
					c	cluster or molecule
					s	supercell
		
		@keyword Atomcount:  atom count (default 0)
		@keyword AtomTypes:  atom types list, orogin, lattice vctors and Geometry array (default None)
		@keyword Origin: 	supercell origin, always ignored but part of .gen file format (default None)
		@keyword Lattice:  lattice vectors array (default None)
		@keyword Geometry: 	atomic positions array (carthesian bohr), for R_0 (default None)
		@keyword AtomLayers: 	list of atom layer index assigments per atom (default None)
		@keyword LayerDict: 	dictionary of (atom layer index: atom layer)-s (default None)
		@keyword AtomCharges: 	list of single atom charges. Autoinitialized to zero if None. (default None)
		@keyword AtomSubTypes:  list of atom subtype strings. Autoinitialized to element symbols if None (default None)
		@keyword LPops:  list of atomic l-shell populations, can be list of lists or None (default None)
		@keyword Axis: inital Dimer axis (default None)
		@keyword noGradInRot: do not perform gradient calculations in rotation step (default False)
		@keyword verbosity: verbosity level, c.f. Verbosity level descriptions in comatsci.constants (default VBL_NORMAL)
		@keyword maxFt: =1e-6 maximum total force for convergence (default 1e-6)
		@keyword maxFp: =1e-6 maximum force in dimer direction for convergence (default 1e-6)
		@keyword maxFtRMS: =1e-8 maximum RMS total force for convergence (default 1e-8)
		@keyword maxFpRMS: =1e-8 maximum RMS force in dimer direction for convergence (default 1e-8)
		@keyword rotskipangle: calculated test angle in radians below which rotation step is skipped (default 0.01)
		@keyword hardTranslation: If!=None, use specified length in Bohr for line search in translation step (default None)
		@keyword maxIt: maximum number of iterations in this run (default 100)
		@keyword fixedAtoms: list of atom indices which are not to be moved (base 0) (default empty list)
		"""
		# check for presence of options dictionary, provide empty one if necessary
		if options==None:
			options={}
		# enforce dimer initialization with Geometry and Axis at the same time
		if ((options.has_key("Geometry") or options.has_key("Axis")) and 
					not (options.has_key("Geometry") and options.has_key("Axis"))):
			raise GeometryError("Dimer midpoint and axis not specified together")
		# unpack options dictionary with defaults
		Mode=options.get("Mode","C")
		Atomcount=options.get("Atomcount",0)
		AtomTypes=options.get("AtomTypes",None) #@UnusedVariable
		Origin=options.get("Origin",None)
		Lattice=options.get("Lattice",None)
		self._R0=options.get("Geometry",None)
		AtomLayers=options.get("AtomLayers",None)
		LayerDict=options.get("LayerDict",None)
		AtomCharges=options.get("AtomCharges",None)
		AtomSubTypes=options.get("AtomSubTypes",None)
		LPops=options.get("LPops",None)
		#Call base class constructor
		Geometry.Geometry.__init__(self, Mode, Atomcount, Origin, Lattice, self._R0, AtomLayers, LayerDict,AtomCharges,AtomSubTypes,LPops)
		# unpack with defaults variables not in base class
		self._DeltaR=options.get("Axis",None)
		self.noGradInRot=options.get("noGradInRot",None)
		self.maxFt=options.get("maxFt",1e-6)
		self.maxFp=options.get("maxFp",1e-6)
		self.maxFtRMS=options.get("maxFtRMS",1e-8)
		self.maxFpRMS=options.get("maxFpRMS",1e-8)
		self.rotSkipAngle=options.get("rotskipangle",0.01)
		self.hardTranslation=options.get("hardTranslation",None)
		self.maxIt=options.get("maxIt",100)
		self.fixedAtoms=options.get("fixedAtoms",[])
		self.verbosity=options.get("verbosity",constants.VBL_NORMAL)
		# set last search direction for CG translation step to None to perform initial steepest descent step
		self.lastFdag=None
		# set counters
		self._setCount=0
		# set atomic constraint defaults
		self._fixedAtoms=[]

#----------------------------------------------------------------------------------------
# fixed atoms property

	def getFixedAtoms(self):
		"""return the list of fixed atoms"""
		return self._fixedAtoms
	
	
	def setFixedAtoms(self,value):
		"""set the list if fixed atom indices from a iterable type
		@param value: the list of atoms to keep fixed"""
		# check if atim indices ar ein range
		for i in value:
			if i >= self.Atomcount:
				raise ValueError("Fixed atom index not in Geometry")
		self._fixedAtoms=value


	fixedAtoms=property(getFixedAtoms,setFixedAtoms,doc="List of atoms to fix during dimer iterations, counts from 0")
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# stepCount read-only property
	@property
	def stepCount(self):
		"""Get number of iterations."""
		return self._setCount
#----------------------------------------------------------------------------------------
	
	
	
	def _reset_derived(self):
		"""reset derived Data, i.e. when reading Dimer info from file
		"""
		#reset energies
		self.E0=None
		self.E1=None
		self.E2=None
		#reset forces
		self.fN=None
		self.f0=None
		self.f1=None
		self.f2=None
		#reset search directions, curvature etc.
		self.Phi=None
		self.C0=None
		self.dC0=None
		#call base class reset method
		Geometry.Geometry._reset_derived(self)
	
	
	
	
	def setR0Geo(self,value):
		"""update center of mass from Geometry Object
		@param value: new center of mass Geometry Object
		"""
		# If R0 has been initialized before, check if geometries are compatible, but ignore origin checking
		# origin is stupid anway
		if self.isInitialized:
			try: 
				self.compatcheck(value)
			except GeometryError,inst:
				if inst.args[0]=='Geometry lattice mismatch':
					if self.verbosity>constants.VBL_SILENCE:
						print("ReactionPath warning: Geometry lattice mismatch")
				else:
					raise
		#now store now c.o.m.
		self.Geometry=copy.deepcopy(value.Geometry)
		self.Mode=copy.deepcopy(value.Mode)
		self.Atomcount=copy.deepcopy(value.Atomcount)
		self.AtomTypes=copy.deepcopy(value.AtomTypes)
		self.Origin=copy.deepcopy(value.Origin)
		self.Lattice=copy.deepcopy(value.Lattice)
		self.AtomLayers=copy.deepcopy(value.AtomLayers)
		self.LayerDict=copy.deepcopy(value.LayerDict)
		self.AtomCharges=copy.deepcopy(value.AtomCharges)
		self.AtomSubTypes=copy.deepcopy(value.AtomSubTypes)
		self.LPops=copy.deepcopy(value.LPops)



	def getR0Geo(self):
		"""return dimer center of mass as Geometry Object
		@return: center-of-mass geometry object
		"""
		if self.isInitialized:
			return self
		else:
			return None
	
	R0=property(getR0Geo,setR0Geo,doc="The dimer center-of-mass geometry object")
	
	
##	# Overwrite base class Geometry attribute
##	
##	def getGeometry(self):
##		""" Return c.o.m. Geometry Array
##		@return: Center of mass coordinartes array"""
##		return self.Geometry
##	
##	
##	
##	def setGeometry(self, value):
##		"""update center of mass from coordinates array
##		@param newCoords: new center of mass coordinates array"""
##		if not shape(self.Geometry)==shape(value):
##			raise GeometryError("Dimer coordinate array mismatch")
##		else:
##			self.Geometry = value
##	Geometry=property(getGeometry,setGeometry,doc="The dimer center-of-mass coordinates array")
	
#----------------------------------------------------------------------------------------
# isInitialized read-only property
	def hasData(self):
		"""check if Dimer configuration is initialized
		@return: true if dimer data is complete
		"""
		if (self.Geometry!=None and self._DeltaR!=None):
			return True
		else:
			return False


	isInitialized=property(fget=hasData,doc="True if dimer data is complete")
#----------------------------------------------------------------------------------------

	
	def setTwoGeometries(self, R1, R2):
		"""set the dimer from two Geometry objects
		@param R1: First geometry
		@param R2: Second geometry
		"""
		# Check if geometries are compatible, but ignore origin checking
		# origin is stupid anway
		try: 
			R1.compatcheck(R2)
		except GeometryError,inst:
				if inst.args[0]=='Geometry lattice mismatch':
					if self.verbosity>constants.VBL_SILENCE:
						print("ReactionPath warning: Geometry lattice mismatch")
				else:
					raise
		# reset the Dimer object
		self._reset_derived()
		# Internally, we store the center-of-mass position and dimer displacement
		# reuse all properties exept coordinate array from R1
		comgeo=(R1.Geometry+R2.Geometry)/2.0
		self.setR0Geo(R1)
		# set proper R0 coordinates and Dimer axis
		self.Geometry=comgeo
		self._DeltaR=((R2.Geometry-R1.Geometry)/2.0).ravel()



##	def setR0DeltaR(self, R0, DeltaR):
##		"""set the dimer from center-of-mass geometry objet and displacement vector
##		@param R0: center-of-mass gemetry
##		@param DeltaR: displacement (R2-R1)/2
##		"""
##		self=R0
##		self.DeltaR=DeltaR



#----------------------------------------------------------------------------------------
# DeltaR property
	def setDeltaR(self,value):
		"""update Delta R
		@param value: the new dimer displacement
		"""
		# check if DeltaR is of compatible dimension to R0 or DeltaR
		if self.Geometry!=None:
			if len(value.ravel())!=len(self.Geometry.ravel()):
				raise ValueError,"Dimension of dimer direction does not fit Geometry"
		elif self._DeltaR!=None:
			if len(value.ravel())!=len(self._DeltaR.ravel()):
				raise ValueError,"Dimension of dimer direction does not existing direction"
		self._DeltaR=value




	def getDeltaR(self):
		"""return DeltaR
		@return: dimer displacement vector
		"""
		if self.isInitialized:
			return self._DeltaR
		else:
			return None


	DeltaR=property(getDeltaR,setDeltaR,doc="The dimer displacement vector")
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# R1 property
	def setR1(self, value):
		"""update R1 of the dimer
		@param value: first dimer geometry
		"""
		# If R0 has been initialized before, check if geometries are compatible, but ignore origin checking
		# origin is stupid anway
		if self.Geometry!=None:
			try: 
				self.compatcheck(value)
			except GeometryError,inst:
				if inst.args[0]=='Geometry lattice mismatch':
					if self.verbosity>constants.VBL_SILENCE:
						print("ReactionPath warning: Geometry lattice mismatch")
				else:
					raise	
		R2=copy.deepcopy(self)
		R2.Geometry-=-self._DeltaR
		self.setTwoGeometries(value,R2)



	def getR1(self):
		"""return R1
		@return: the first dimer geometry
		"""
		if self.isInitialized:
			R1=copy.deepcopy(self)
			R1.setcoordinates(self.Geometry+numpy.reshape(self._DeltaR,numpy.shape(self.Geometry)))
			return R1
		else:
			return None


	R1=property(getR1,setR1,doc="The first geometry of the dimer")
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# R2 property
	def setR2(self, value):
		"""update R2 of the dimer
		@param value: first dimer geometry
		"""
		# If R0 has been initialized before, check if geometries are compatible, but ignore origin checking
		# origin is stupid anway
		if self.Geometry!=None:
			try: 
				self.Geometry.compatcheck(value)
			except GeometryError,inst:
				if inst.args[0]=='Geometry lattice mismatch':
					if self.verbosity>constants.VBL_SILENCE:
						print("ReactionPath warning: Geometry lattice mismatch")
				else:
					raise
		R1=copy.deepcopy(self.Geometry)
		R1.Geometry+=numpy.reshape(self._DeltaR,numpy.shape(R1.Geometry))
		self.setTwoGeometries(R1,value)



	def getR2(self):
		"""return R1
		@return: the first dimer geometry
		"""
		if self.isInitialized:
			R2=copy.deepcopy(self)
			R2.setcoordinates(self.Geometry-numpy.reshape(self._DeltaR,numpy.shape(self.Geometry)))
			return R2
		else:
			return None


	R2=property(getR2,setR2,doc="The first geometry of the dimer")
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# direction property
	def getDimerDir(self):
		"""return normalized dimer displacement vector
		@return: normalized DeltaR
		"""
		if self.isInitialized:
			return self._DeltaR/math.sqrt(numpy.dot(self._DeltaR.ravel(),self._DeltaR.ravel()))
		else:
			return None




	def setDimerDir(self,value):
		"""set direction of displacement vector to direction Direction, keeping displacement length
		@param value: vector defining new displacement direction
		"""
		#first, check if we are initialized, otherwise this operation makes no sense
		if not self.isInitialized:
			raise "Trying to set direction on an uninitialized dimer"
		#now check if Direction vector is of suitable dimension
		if len(value.ravel())!=len(self._DeltaR.ravel()):
			raise ValueError,"Dimer direction vector incompatible to existing dimer data"
		#if we came this far, store displacement length, normalize Direction and store new displacement vector
		length=self.displacement
		value/=math.sqrt(numpy.dot(value,value))
		self._DeltaR=value*length


	direction=property(getDimerDir,setDimerDir,doc="The dimer displacement direction")
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# displacement property
	def getDisplacement(self):
		"""return dimer displacement length
		@return: length of DeltaR
		"""
		if self.isInitialized:
			return math.sqrt(numpy.dot(self._DeltaR.ravel(),self._DeltaR.ravel()))
		else:
			return None




	def setDisplacement(self,value):
		"""set dimer displacement length, keeping displacement direction
		@param value: new dimer displacement vector length
		"""
		#check if dimer is initialized, operation makes no sense otherwise
		if not self.isInitialized:
			raise "Trying to set displacement on an uninitialized dimer"
		#scale DeltaR to new length
		self._DeltaR*=float(value)/math.sqrt(numpy.dot(self._DeltaR.ravel(),self._DeltaR.ravel()))


	displacement=property(getDisplacement,setDisplacement,doc="The dimer displacement vector length")
#----------------------------------------------------------------------------------------


	def rotate(self, rotdir, angle):
		"""rotate the by angle radians into direction. <b>direction is not the axis of rotation, together with the dimer axis, it spans the plane of rotation</b>.
		@param rotdir: initial tangent vector of rotation. Is internally orthoganalized to dimer axis and normalized
		@param angle: angle of rotation in radians
		"""
		#check if dimer is initialized, operation makes no sense otherwise
		if not self.isInitialized:
			raise "Trying to rotate an uninitialized dimer"
		#check, if rotdir is of equal dimension to dimer
		if not len(self.direction.ravel())==len(rotdir.ravel()):
			raise ValueError,"Tring to rotate dimer into direction of unequal dimension"
		#Gram-Schmidt orthoganlize rotdir to the dimer axis direction
		#remember that self.direction is normalized
		rotdir-=(numpy.dot(self.direction,rotdir)*self.direction)
		#normalize rotdir
		rotdir/=math.sqrt(numpy.dot(rotdir,rotdir))
		# now rotate
		d0=self.direction
		# d'=d0*cos(phi)+rotdir*sin(phi)
		self.setDimerDir(d0*math.cos(angle)+rotdir*math.sin(angle))
##		#prepare some intermediate values to make formula more readable
##		Rsq=self.displacement**2
##		Rpwsix=Rsq**3
##		cosinephi=math.cos(angle)
##		#now rotate dimer axis
##		self._DeltaR+=math.sqrt(abs(((Rpwsix)/cosinephi**2)-Rsq))*rotdir
##		self._DeltaR*=(cosinephi/Rsq)



	def calcEnergiesForces(self, scheduler, charge=0):
		"""Calculate energies and forces at R0 and R1 using primed scheduler, extrapolate f2 and E2 after J. Phys. Chem. 123 224101
		This implementation will do one unneccesary E1/f1 calculation after convergence , as is is checked only after E0/f0 and E1/f1 have been calculated
		@param scheduler: scheduler used to distribute the two jobs
		@param charge: total charge of the system in electrons (default=0)
		"""
		if not self.isInitialized:
			raise GeometryError("Trying to calculate Energy and force on uninitialized Dimer")
		else:
			# assemble options dictionaries
			joblist=[]
			joblist.append({
				"Geometry":copy.deepcopy(self.R0),
				"Charge":charge,
				"steplabel":"R0",
			})
			joblist.append({
				"Geometry":copy.deepcopy(self.R1),
				"Charge":charge,
				"steplabel":"R1",
			})
			# let the scheduler perform the jobs on the list
			results=scheduler.perform(joblist)
			# now unpack the results
			self.E0=results[0]["Energy"]
			self.f0=self.atomPointConstraint(results[0]["Forces"].ravel())
			self.E1=results[1]["Energy"]
			self.f1=self.atomPointConstraint(results[1]["Forces"].ravel())
			# Calculate the necessary properties prior to rotation
			self.__extrapolateEfC()




	def __extrapolateEfC(self):
		"""Perform the calculations necessary after calling the external calculator and before the rotatation step (second box from top in Fig. 2 of J. Chem. Phys. 123, 224101)
		"""
		# extrapolate f2
		self.f2=2.0*(self.f0-self.f1)
		# calculate constrained rotational force and normalize to dimer axis
		self.fN=self.f1-self.f2
		# Gram-Schmidt orthoganlize fN to the dimer axis direction
		# (remember that self.direction is normalized)
		self.fN-=(numpy.dot(self.direction,self.fN)*self.direction)
		# rotational search direction is normalized rotational force vector
		self.Phi=self.fN/math.sqrt(numpy.dot(self.fN,self.fN))
		# extrapolate E2 using proper formula, depending on rotation gradient calculation setting
		if self.noGradInRot:
			self.E2=2.*self.E1-self.E0+3.*self.displacement*numpy.dot(self.direction,self.f0)
			self.C0=2.*(self.E1 -self.E0+numpy.dot(self.f0,self.direction)*self.displacement)/self.displacement**2
		else:
			self.E2=2.*(self.E0-0.25*self.displacement*numpy.dot(self.f1-self.f2,self.direction))-self.E1
			self.C0=0.5*numpy.dot(self.f2-self.f1,self.direction)/2.*self.displacement
		self.dC0=-numpy.dot(self.fN,self.Phi)/self.displacement
	
	
	
	
	def _rotationStep(self,scheduler):
		"""rotate dimer into lowest curvature direction"""
		#save old position
		self.oldDeltaR=self.DeltaR
		#
		if self.noGradInRot:
			tempDimer=copy.deepcopy(self)
			rotangle=constants.PI/4.
			tempDimer.rotate(self.Phi,rotangle)
			Er,fr=tempDimer.calcE1(scheduler)
			Cr=2*(Er-self.E0+numpy.dot(self.DeltaR,tempDimer.direction))/self.displacement**2
		else:
			rotangle=0.5*math.atan(-self.dC0/(2.*self.C0))
			# if calculated phi_1 is smaller than a predefined threshold, skip the rotation step
			if rotangle<self.rotSkipAngle:
				return
			tempDimer=copy.deepcopy(self)
			tempDimer.rotate(self.Phi,rotangle)
			Er,fr=tempDimer.calcE1(scheduler)
			Cr=0.5*numpy.dot(self.f2-fr,self.direction)/self.displacement
		# determine fourier coefficients of E(phi) and phi_min
		a1=(self.C0-Cr+0.5*self.dC0*math.sin(2.*rotangle))/(1-math.cos(2.*rotangle))
		b1=0.5*self.dC0
##		a0=2.*(self.C0-a1) # not needed
		phi_min=0.5*math.atan(b1/a1)
		# rotate self into minimum curvature mode
		self.rotate(self.Phi,phi_min)




	def _translationStep(self, scheduler):
		"""translation step of dimer algorithm"""
		#determine modified force at original position
		if self.C0>0:
			fdag_0=-(numpy.dot(self.f0,self.direction)*self.direction)
		else:
			fdag_0=self.f0-2.*(numpy.dot(self.f0,self.direction)*self.direction)
		# *** determine search direction ***
		# start with simple SD
		G=fdag_0
		# if last direction is stored, calculate conjugate search direction using
		# Polak-Ribiere formula
		if self.lastFdag!=None:
			G-=self.lastFdag*(numpy.dot(fdag_0,fdag_0-self.lastFdag)/numpy.dot(self.lastFdag,self.lastFdag))
		# normalize translation direction
		G/=math.sqrt(numpy.dot(G,G))
		# *** end search direction calculation***
		#determine trial translation length, calculate if no hard default is provided
		if self.hardTranslation==None:
			translation=0.5*numpy.dot(G,fdag_0)/abs(self.C0)
		else:
			translation=self.hardTranslation
		# translate trial dimer and calculate forces
		testDimer=copy.deepcopy(self)
		testDimer.Geometry+=numpy.reshape((G*translation),numpy.shape(testDimer.Geometry))
		Nrg, f_new = testDimer.calcE0(scheduler) #@UnusedVariable
		# calculate modified force at trial position
		if self.C0>0:
			fdag_t=-(numpy.dot(self.f0,self.direction)*self.direction)
		else:
			fdag_t=testDimer.f0-2.*(numpy.dot(testDimer.f0,self.direction)*self.direction)
		# *** determine root of assumed linear force in search direction ***
		# f(r0)=a*r0+f0==0
		# project modified forces onto search direction
		fs0=numpy.dot(fdag_0,G)
		fst=numpy.dot(fdag_t,G)
		# determine offset and slope of linear force function
		# {a=(f1-f2)/(r1-r2)
		a=(fs0-fst)/(-translation)
		# {fo=f2-a*r2}
		f0=fst-a*translation
		# {r0=-f0/a}
		r0=-f0/a
		# *** end linear root search ***
		# translate dimer to assumed root position
		self.Geometry+=numpy.reshape((G*r0),numpy.shape(self.Geometry))
		# all previously calculated data is now invalid, so discard it!
		self._reset_derived()
		# translation step finished




	def calcE1(self, scheduler, charge=0):
		"""calculate E1 and f1 for dimer rotation step
		@return: (E1,f1) tuple of energy and force at R1"""
		joblist=[]
		joblist.append({
			"Geometry":copy.deepcopy(self.R1),
			"Charge":charge,
			"steplabel":"R1",
		})
		# let the scheduler perform the jobs on the list
		results=scheduler.perform(joblist)
		# return results
		return (results[0]["Energy"],self.atomPointConstraint(results[0]["Forces"].ravel()))
	
	
	
	
	
	def calcE0(self,scheduler,charge=0):
		"""calculate E0 and f0 for dimer rotation step
		@return: (E0,f0) tuple of energy and force at R1"""
		joblist=[]
		joblist.append({
			"Geometry":copy.deepcopy(self.R0),
			"Charge":charge,
			"steplabel":"R0",
		})
		# let the scheduler perform the jobs on the list
		results=scheduler.perform(joblist)
		# return results
		return (results[0]["Energy"],self.atomPointConstraint(results[0]["Forces"].ravel()))
	
	
	
	
	
	def readfmg(self, filename):
		"""read dimer from fmg file, using minidom
		@param filename: input file name
		"""
		# Read the input file into a string
		inFile=utils.compressedopen(filename,"r")
		lines=list(inFile)
		xmlString="".join(lines)
		inFile.close()
		del lines
##		dom = xml.dom.minidom.parseString(DIMER_DTD+xmlString)
		dom = xml.dom.minidom.parseString(xmlString)
		self.handleDimerDom(dom.getElementsByTagName("dimer")[0])
		#clear up memory after reading
		dom.unlink()
	
	
	
	
	def handleDimerDom(self,dimerDom):
		"""handle dimer DOM constructed by input file parser
		@param dimerDom: the input file DOM object to handleNameClass
		"""
		#Constant tables for unit attributes
		lunits={"ang":constants.ANGSTROM, "au":constants.BOHR}
		eunits={"H":constants.HARTREE, "au":constants.HARTREE, "eV":constants.EVOLT}
		#clear derived info and set some default values
		self._reset_derived()
		self.NoGradInRot=False
		#handle geometry info using base class handler
		ingeo=dimerDom.getElementsByTagName("geometry")[0]
		self.handlegeometry_dom(ingeo)
		# *** handle dimer Separation vector ***
		DeltaR=dimerDom.getElementsByTagName("DeltaR")[0]
		if DeltaR.hasAttribute("lunit"):
				lunit=lunits[DeltaR.getAttribute("lunit").strip()]
		else:
			#default length unit is Angstrom
			lunit=lunits["ang"]
		temp=[ float(t)/lunit for t in (DeltaR.childNodes[0].data.strip().split()) ]
		self._DeltaR=numpy.array(temp,dtype=float)
		# *** end dimer separation vecot handling ***
		#get Gradient Calc in Rotation Step boolean (False is already set as default)
		if len(dimerDom.getElementsByTagName("NoGradInRot"))>0:
			self.noGradInRot=True
		# *** Energies and forces handling ***
		#if a complete set is present, handle E0,E1,f0,f1, otherwise skip all calculated Data
		E0=dimerDom.getElementsByTagName("E0")
		E1=dimerDom.getElementsByTagName("E1")
		f0=dimerDom.getElementsByTagName("f0")
		f1=dimerDom.getElementsByTagName("f1")
		if len(E0)>0 and len(E1)>0 and len(f0)>0 and len(f1)>0:
			if E0[0].hasAttribute("eunit"):
				eunit=eunits[E0[0].getAttribute("eunit").strip()]
			else:
				#default energy unit is Hartrees
				eunit=eunits["au"]
			self.E0=float(E0[0].childNodes[0].data.strip())/eunit
			if E1[0].hasAttribute("eunit"):
				eunit=eunits[E1[0].getAttribute("eunit").strip()]
			else:
				#default energy unit is Hartrees
				eunit=eunits["au"]
			self.E1=float(E1[0].childNodes[0].data.strip())/eunit
			# we only support H/Bohr as force unit
			self.f0=numpy.array([ float(tmp) for tmp in (f0[0].childNodes[0].data.strip().split()) ] )
			self.f1=numpy.array([ float(tmp) for tmp in (f1[0].childNodes[0].data.strip().split()) ] )
			#E2,f2,fN and all curvature data from intput file is always ignored and re-calculated instead!
			self.__extrapolateEfC()
		# *** end Energies and forces handling ***




	def checkConvergence(self):
		"""Check converence of dimer saddle-point search
		@return: boolean overall convergence result"""
		#initialize return value, if any test fails, report non-convergence
		converged=True
		#calculate f0 in dimer direction
		f0_proj=self.f0*(numpy.dot(self.f0,self.direction)/math.sqrt(numpy.dot(self.f0,self.f0)))
		#check maximum force
		if max(self.f0.ravel())>self.maxFt:
			converged=False
			if self.verbosity>=constants.VBL_NORMAL:
				print("max f0 component            : No")
		else:
			if self.verbosity>=constants.VBL_NORMAL:
				print("max f0 component            : Yes")
		#check maximum force in dimer direction
		if max(f0_proj.ravel())>self.maxFp:
			converged=False
			if self.verbosity>=constants.VBL_NORMAL:
				print("max projected f0 component  : No")
		else:
			if self.verbosity>=constants.VBL_NORMAL:
				print("max projected f0 component  : Yes")
		#check RMSD force
		if math.sqrt(numpy.dot(self.f0,self.f0)/float(len(self.f0)))>self.maxFtRMS:
			converged=False
			if self.verbosity>=constants.VBL_NORMAL:
				print("RMS f0                      : No")
		else:
			if self.verbosity>=constants.VBL_NORMAL:
				print("RMS f0                      : Yes")
		#check projected RMSD force
		if math.sqrt(numpy.dot(f0_proj,f0_proj)/float(len(self.f0)))>self.maxFpRMS:
			converged=False
			if self.verbosity>=constants.VBL_NORMAL:
				print("RMS projected f0            : No")
		else:
			if self.verbosity>=constants.VBL_NORMAL:
				print("RMS projected f0            : Yes")
		return converged




	def readfile(self, filename,typespec=None):
			"""read geometry from specified filename
			@param filename: input file name
			@param typespec: lowercase string, specifying file type, autodetect if omitted
			known types:
			<ul>
			<li>fmg</li>
			</ul>
			"""
			typefunctions={
				'fmg': self.readfmg
				}
			if typespec!=None:
				ftype=typespec.strip('.').strip().lower()
			else:
				ftype=filename.strip()[-3:].lower()
			if ftype in typefunctions.keys():
				typefunctions[ftype](filename)
			else:
				raise GeometryError("Unknown dimer input file type")



##	def _buildOutputDOM(self):
##		"""builds xml DOM object for output
##		@return: DOM representation for .fmg output"""
##		dom=xml.dom.minidom



	def getFmgString(self):
		"""build fmg string representation of the dimer object
		@return: string containing .fmg format dimer declaration"""
		#start with the dimer root element
		fmgstring="<dimer>\n"
		#DTD dmeands geometry to be the first element in dimer
		fmgstring+=Geometry.Geometry.getFmgString(self)+"\n"
		#Delta_R is next
		fmgstring += self._array_element(self.DeltaR,"DeltaR")
		#noGradInRot boolean
		if self.noGradInRot:
			fmgstring+="<NoGradInRot />\n"
		#Energies, if present
		if self.E0!=None and self.E1!=None:
			fmgstring+="<E0>{0:18.12e}</E0>\n".format(self.E0)
			fmgstring+="<E1>{0:18.12e}</E1>\n".format(self.E1)
			fmgstring+="<E2>{0:18.12e}</E2>\n".format(self.E2)
		#forces if present
		if self.f0!=None and self.f1!=None:
			fmgstring += self._array_element(self.f0,"f0")
			fmgstring += self._array_element(self.f1,"f1")
			fmgstring += self._array_element(self.f2,"f2")
			fmgstring += self._array_element(self.fN,"fN")
		#curvature if present
		if self.C0!=None:
			fmgstring+="<curvature>{0:18.12e}</curvature>\n".format(self.C0)
		#close dimer element and return
		return fmgstring+"</dimer>"




	def _array_element(self, data, tag, indentlevel=0):
		"""return a string containing an xml element of type tag
		containing data as a (-1,3) Block of floats
		@param data: the array to be included in the element
		@param tag: name of the xml element to generate
		@indentlevel=0 number of tab-stops to insert in front of all lines"""
		returnString=("\t"*indentlevel)+("<{0:s}>\n".format(tag))
		tmp=numpy.reshape(data,(-1,3))
		for i in tmp:
			returnString+=("\t"*(indentlevel+1))+("{0:18.12e}\t{1:18.12e}\t{2:18.12e}\n".format(i[0],i[1],i[2]))
		returnString+=("\t"*indentlevel)+("</{0:s}>\n".format(tag))
		return returnString

	
	fmgString=property(getFmgString,doc="String representation of dimer in .fmg (xml) format")



	def dimerIterate(self, scheduler):
		"""iterate improved dimer method
		@param scheduler: scheduler instance to use for calculations
		"""
		# some general chatter
		if self.verbosity>=constants.VBL_NORMAL:
			print("""
     Number of atoms        : {0:d}
     Number of fixed atoms  : {1:d}
     Dimer separation       : {2:e} Ang
			""".format(self.Atomcount,len(self.fixedAtoms),self.displacement*constants.ANGSTROM))
			print("initial energies and forces calculation")
		# first energies and forces calculation (if not read from file)
		if self.E0==None:
			self.calcEnergiesForces(scheduler)
			# write convergence information into a new file
			dimerInfoFile=open("dimerinfo.dat","w")
			print("#iteration E0 max(f_0) RMS (f_0) C0",file=dimerInfoFile)
			dimerInfoFile.flush()
		else:
			if self.verbosity>=constants.VBL_NORMAL:
				print("Using Energies and forces from previous run")
			# append convergence information to existing file
			dimerInfoFile=open("dimerinfo.dat","a")
		self.writexyz("dimer-animation.xyz","w")
		if self.verbosity>=constants.VBL_NORMAL:
			print("E0                          = {0:12.6f} H".format(self.E0))
			print("max f0                      = {0:12.6f} H/Bohr".format(max(self.f0.ravel())))
			print("RMS f0                      = {0:12.6f} H/Bohr".format(math.sqrt(numpy.dot(self.f0,self.f0)/float(len(self.f0)))))
			print("C0                          = {0:12.6f} H/Bohr**2".format(self.C0))
		sys.stdout.flush()
		# write initial data to convergence info file
		print("{0:d}\t{1:12.8e}\t{2:12.8e}\t{3:12.8e}\t{4:12.8e}".format(0,self.E0,max(self.f0.ravel()),math.sqrt(numpy.dot(self.f0,self.f0)/float(len(self.f0))),self.C0),file=dimerInfoFile)
		dimerInfoFile.flush()
		# check if we can skip iterations (you wish)
		if self.checkConvergence():
			print("Already converged, nothing to optimize.")
			return
		if self.verbosity>=constants.VBL_NORMAL:
			print("starting improved dimer iterations")
			sys.stdout.flush()
		for i in range(self.maxIt):
			if self.verbosity>=constants.VBL_NORMAL:
				print("iteration {0:5d}".format(i+1))
				print("rotation step")
				sys.stdout.flush()
			self._rotationStep(scheduler)
			if self.verbosity>=constants.VBL_NORMAL:
				print("translation step")
				sys.stdout.flush()
			self._translationStep(scheduler)
			self.writefmg("dimercheckpoint.fmg")
			if self.verbosity>=constants.VBL_NORMAL:
				print("recalculating energies and forces")
				sys.stdout.flush()
			self.calcEnergiesForces(scheduler)
			self.writefmg("dimercheckpoint.fmg")
			self._setCount+=1
			# write data to convergence info file
			print("{0:d}\t{1:12.8e}\t{2:12.8e}\t{3:12.8e}\t{4:12.8e}".format(i+1,self.E0,max(self.f0.ravel()),math.sqrt(numpy.dot(self.f0,self.f0)/float(len(self.f0))),self.C0),file=dimerInfoFile)
			dimerInfoFile.flush()
			self.writexyz("dimer-animation.xyz","a")
			if self.verbosity>=constants.VBL_NORMAL:
				print("E0                          = {0:12.6f} H".format(self.E0))
				print("max f0                      = {0:12.6f} H/Bohr".format(max(self.f0.ravel())))
				print("RMS f0                      = {0:12.6f} H/Bohr".format(math.sqrt(numpy.dot(self.f0,self.f0)/float(len(self.f0)))))
				print("C0                          = {0:12.6f} H/Bohr**2".format(self.C0))
				sys.stdout.flush()
				print("convergence check")
			# check for convergence and stop if converged
			if self.checkConvergence():
				dimerInfoFile.close()
				return
			# check for stop file and stop if present
			if os.path.exists("STOP_PESTO"):
				if self.verbosity>=constants.VBL_QUIET:
					print("Stop-file encountered, stopping dimer iterations.")
				dimerInfoFile.close()
				return
		if self.verbosity>=constants.VBL_QUIET:
			dimerInfoFile.close()
			print("maximum number of iterations reached. stopping.")
			sys.stdout.flush()



	def atomPointConstraint(self, inforce):
		"""fix the specified atom's positions
		@param inforce: the force vector to apply contraints to
		@return: constrained forces"""
		# now shape to (-1,3) and zero out the specified atomic forces
		outforce=numpy.reshape(inforce,(-1,3))
		zeroforce=numpy.array([0.,0.,0.])
		for i in self.fixedAtoms:
			outforce[i]=zeroforce
		# put output array in same shape as input and return
		return numpy.reshape(outforce,numpy.shape(inforce))
