## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# NEBPath.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
from __future__ import print_function
import numpy

import os, sys, copy

from .. import geometry,constants

from .reactionpath import Reactionpath

#TODO: for 0.5.0 Add ABNR NEB optimizer

class NEBPath(Reactionpath):
	"""Represents a reaction Path for NEB relaxation"""

	def __init__(self,icheckpointdir,ifixedatoms,
		istepwidth,iforcemode,itangmode,iclimber,ispringk,
		icmode,ifmax,ifrms,imaxit,irelmode,charge=0,verbosity=None):
		"""initialise NEBPath
		@param icheckpointdir: c.f. Reactionpath documentation
		@param ifixedatoms: c.f. Reactionpath documentation
		@param istepwidth: initial stepwidth for velocity verlet
		@param iforcemode: NEB force calculation mode
			* B{s} standard NEB
			* B{c} climbing image NEB
		@param itangmode: NEB tangent calculation mode
			* B{s} standard tangents
			* B{w} weighted tangents after JPC 113,9978(2000)
			* B{p} spline tangents
		@param iclimber: climbing image for climbing image NEB
		@param ispringk: NEB spring constant
		@param icmode: c.f. L{Reactionpath} documentation
		@param ifmax: c.f. L{Reactionpath} documentation
		@param ifrms: c.f. L{Reactionpath} documentation
		@param imaxit: c.f. L{Reactionpath} documentation
		@param irelmode: relaxation mode (ignored for now)<ul>
			* B{v} projected velocity verlet
			* B{s} adaptive displacement steepest descent
		@param charge: system charge
		@param verbosity: Verbosity level. Choose NORMAL as default. (default None)
		"""
		#base class constructor
		Reactionpath.__init__(self,icheckpointdir,ifixedatoms,icmode,ifmax,ifrms,imaxit,charge,verbosity)
		if self.verbosity>=constants.VBL_NORMAL:
			print("Initializing NEB Path.")
		self.nebforces=[]
		self.tangents=[]
		self.velocities=[]
		self.amasses=[]
		if iforcemode=='s':
			self.forcefunc=self.stdnebforces
			if self.verbosity>=constants.VBL_NORMAL:
				print("Will use standard NEB mode")
			self.convcheckfunc=self._sconvcheck
		elif iforcemode=='c':
			if iclimber==-1:
				raise("invalid climber specified for climbing image NEB")
			else:
				self.climber=iclimber
			self.forcefunc=self.cinebforces
			self.convcheckfunc=self._ciconvcheck
			if self.verbosity>=constants.VBL_NORMAL:
				print("Will use climbing image NEB mode")
				print("Image No. {0:d} selected as climbing image".format(self.climber))
		else:
			raise ValueError("Unknown NEB force mode")
		if itangmode=='w':
			self.tangfunc=self._calcweightedtangents
			if self.verbosity>=constants.VBL_NORMAL:
				print("Will use weighted tangents")
		elif itangmode=='s':
			self.tangfunc=self._calctangents
			if self.verbosity>=constants.VBL_NORMAL:
				print("Will use simple tangents")
		elif itangmode=='p':
			self.tangfunc=self._calcSplineTangents
			if self.verbosity>=constants.VBL_NORMAL:
				print("Will use spline tangents")
		else:
			raise ValueError("Unknown NEB tangents mode")
		self.springk=ispringk
		if self.verbosity>=constants.VBL_NORMAL:
			print("Spring constant set to {0:f}".format(self.springk))
		if irelmode=="v":
			self.relaxor=self.veloverlet
			if self.verbosity>=constants.VBL_NORMAL:
				print("Relaxing path using projected velocity verlet")
			self.minstepsize=0.0001
		elif irelmode=="s":
			self.relaxor=self.adsd
			self.minstepsize=0.00001
			if self.verbosity>=constants.VBL_NORMAL:
				print("Relaxing path using adaptive displacement steepest descent")
		else:
			raise ValueError("Unknown NEB relaxation mode")
		self.dt=istepwidth
		if self.verbosity>=constants.VBL_NORMAL:
			print("Initial stepwidth set to {0:6.3f}".format(self.dt))
		self.dtadapted = 0 # in case of overshoot, no further stepsize adaption should take place
		self.__ovr_history_available=False
		self.__forces_history_available=False
		self.__ovr_intermediate_available=False
		# path debugging mode, to find problems with the reactionopath being relaxed
		self._writePathDebugInfo=False
		self._writePathDebugInfoEnergiesFile=None

	
	
	#writePathDebugInfo methods and property
	def setwritePathDebugInfo(self, value):
		self._writePathDebugInfo=bool(value)
	
	def getwritePathDebugInfo(self):
		return self._writePathDebugInfo
	
	writePathDebugInfo=property(fget=getwritePathDebugInfo, fset=setwritePathDebugInfo,fdel=None,doc="If true, write path debug information in every NEB step")


	
	def _calcSplineTangents(self):
		"""calculate local tangents at all path points using the spline representation of the path
		"""
		# always re-initialize spline representation
		self._genSplineRep()
		#now generate tangents for all images from coordinates vector spline
		for i in range(self.numImages):
			position=self.splineRep["geo"]._nodesx[i]
			temp=self.splineRep["geo"].splder(position)
			tau=temp/numpy.dot(temp,temp)
			self.tangents.append(tau)
	
	

	def _calctangents(self):
		"""Calculate local tangents at all path points, returns 0 tangents for start and endpoint"""
		self.tangents=[]
		self.tangents.append(numpy.zeros((3*self.Atomcount),dtype=float))
		for i in range(1,(self.numimages()-1)):
			Rminus=self.geos[i-1].Geometry.ravel()
			R=self.geos[i].Geometry.ravel()
			Rplus=self.geos[i+1].Geometry.ravel()
			tauleft=R-Rminus
			tauleft/=numpy.sqrt(numpy.dot(tauleft,tauleft))
			tauright=Rplus-R
			tauright/=numpy.sqrt(numpy.dot(tauright,tauright))
			tau=tauleft+tauright
			temp=numpy.sqrt(numpy.dot(tau,tau))
			if temp!=0:
				tau/=temp
			self.tangents.append(tau)
		self.tangents.append(numpy.zeros((3*self.Atomcount),dtype=float))



	def _calcweightedtangents(self):
		"""Calculate Energy-weighted local tangents Following JPC 113,9978(2000), 
		if no energies are defined, do standard tangent calculation"""
		if self.energies == []:
			self.calctangents()
		else:
			self.tangents=[]
			self.tangents.append(numpy.zeros((3*self.Atomcount),dtype=float))
			for i in range(1,(self.numimages()-1)):
				Rminus=self.geos[i-1].Geometry.ravel()
				R=self.geos[i].Geometry.ravel()
				Rplus=self.geos[i].Geometry.ravel()
				tauleft=R-Rminus
				tauright=Rplus-R
				if self.energies[i+1] > self.energies[i] and self.energies[i] > self.energies[i-1]:
					tau=tauright
				elif self.energies[i+1] < self.energies[i] and self.energies[i] < self.energies[i-1]:
					tau=tauleft
				elif self.energies[i+1] > self.energies[i-1]:
					tau=(tauright*numpy.maximum(abs(self.energies[i+1]-self.energies[i]),abs(self.energies[i-1]-self.energies[i])))
					tau+=(tauleft*numpy.minimum(abs(self.energies[i+1]-self.energies[i]),abs(self.energies[i-1]-self.energies[i])))
				else:
					tau=(tauleft*numpy.maximum(abs(self.energies[i+1]-self.energies[i]),abs(self.energies[i-1]-self.energies[i])))
					tau+=(tauright*numpy.minimum(abs(self.energies[i+1]-self.energies[i]),abs(self.energies[i-1]-self.energies[i])))
				temp=numpy.sqrt(numpy.dot(tau,tau))
				if temp!=0:
					tau/=temp
				self.tangents.append(tau)
			self.tangents.append(numpy.zeros((3*self.Atomcount),dtype=float))



	def _remtangentforces(self,forces,tangents):
		"""Remove Tangential compontent from forces, used in NEB-Froces methods
		True Forces must have been calculated beforehand!
		@param forces: forces array to remove tangent components frombuffer
		@param tangents: array storing tangent directions (must be normalized to unity)
		"""
		tmp=[]
		for i in range(self.numimages()):
			absforce=numpy.sqrt(numpy.dot(forces[i].ravel(),forces[i].ravel()))
			forcedir=forces[i].ravel()/absforce
			fdottau=numpy.dot(forcedir.ravel(),tangents[i])
			tmp.append(numpy.subtract(forces[i],numpy.reshape(((fdottau*absforce)*tangents[i]),numpy.shape(forces[i]))))
		return tmp



	def realnormalforces(self):
		"""@return: the real forces normal to the path direction"""
		rnf = self._remtangentforces(self.realforces,self.tangents)
		return rnf



	def realtangentforces(self):
		"""@return: the real forces parallel to the path direction"""
		returnfrc=[]
		for i in range(self.numimages()):
			absforce=numpy.sqrt(numpy.dot(self.realforces[i].ravel(),self.realforces[i].ravel()))
			forcedir=self.realforces[i].ravel()/absforce
			fdottau=numpy.dot(forcedir.ravel(),self.tangents[i])
			returnfrc.append(self.tangents[i]*fdottau*absforce)
		return returnfrc



	def stdnebforces(self):
		"""Return NEB Forces
		True forces must have been calculated beforehand!"""
		# first remove tangential components from true force
		self.tangfunc()
		tempforces=self._remtangentforces(self.realforces,self.tangents)
		#now calculate spring forces
		for i in range(1,self.numimages()-1):
			Rminus=self.geos[i-1].Geometry.ravel()
			R=self.geos[i].Geometry.ravel()
			Rplus=self.geos[i+1].Geometry.ravel()
			tauleft=R-Rminus
			tauright=Rplus-R
			springforce=self.springk*(numpy.sqrt(numpy.dot(tauright,tauright))-numpy.sqrt(numpy.dot(tauleft,tauleft)))
			numpy.add(tempforces[i],springforce*numpy.reshape(self.tangents[i],numpy.shape(tempforces[i])),tempforces[i])
		return tempforces



	def cinebforces(self):
		"""Return Climbing image NEB forces.
		True forces must have been calculated beforehand!
		calls stdnebforces"""
		tempforces = self.stdnebforces()
		forceshape = numpy.shape(tempforces[self.climber])
		tangentforce = self.realtangentforces()[self.climber].ravel()
		normalforce = self.realnormalforces()[self.climber].ravel()
		climbforce = normalforce-tangentforce
		climbforce = numpy.reshape(climbforce,forceshape)
		tempforces[self.climber] = climbforce
		return tempforces



	def initveloverlet(self):
		"""initialize velocity-verlet relaxation scheme"""
		for i in range(self.numimages()): #@UnusedVariable
			self.velocities.append(numpy.zeros((3*self.Atomcount),dtype=float))



	def veloverlet(self,forces):
		"""Do one velocity-verlet step
		@param forces: forces vector
		"""
		for i in range(1,self.numimages()-1):
			forcedir=forces[i].ravel()/numpy.sqrt(numpy.dot(forces[i].ravel(),forces[i].ravel()))
			if numpy.dot(self.velocities[i].ravel(),forcedir) >= 0:
				self.velocities[i]=forcedir*(numpy.dot(self.velocities[i],forcedir))
			else:
				self.velocities[i]=numpy.zeros((3*self.Atomcount),dtype=float)
				if self.verbosity>=constants.VBL_QUIET:
					print("Warning: Verlet velocity vector for image {0:d} opposed to force.".format(i))
			accel=0.5*forces[i]
			for j in range(self.Atomcount):
				accel[j]/=self.amasses[j]
			self.velocities[i]+=accel.ravel()*self.dt
			coords=self.geos[i].Geometry
			displace=(numpy.reshape(self.velocities[i],numpy.shape(accel))+accel)*self.dt
			for j in range(self.Atomcount):
				d=numpy.sqrt(numpy.dot(displace[j],displace[j]))
				if d > 0.4:
					if self.verbosity>=constants.VBL_QUIET:
						print("Warning: Displacement > 0.4 Ang detected in image {0:d}.".format(i))
			coords+=displace
			self.geos[i].setcoordinates(coords)
		# kill possibly stored spline representation
		self.splineRep=None


	def _abnrobjective(self):
		"""This is the objective function for abnr relaxation
		after JCP _119_ 12708 (2003)"""
		# First half of S: \sum_{i=0}^{N} V(r_i)
		objective=0
		for i in self.energies:
			objective+=i
		# second half of S: \sum_{i=1}^{N-1} 1/2 k (\Deltal_i^{RMS}-\Deltal^{RMS})^2
		# first calculate \Deltal^{RMS}
		delta=0
		for i in range(1,self.numimages()-1):
			delta+=self.geos[i].Geometry.ravel()-self.geos[i-1].Geometry.ravel()
		delta/=(self.numimages()*self.Atomcount)
		DeltalRMS=numpy.sqrt(delta)
		# now calculate \Deltal_i^{RMS}
		Deltal_i=[0] # i runs from 1, insert a dummy i=0 here
		for i in range(1,self.numimages()-1):
			delta=self.geos[i].Geometry.ravel()-self.geos[i-1].Geometry.ravel()
			Deltal_i.append(numpy.sqrt(delta/self.Atomcount))
		# add the second sum to out objective
		for i in range(1,self.numimages()-1):
			objective+=0.5*self.springk*((Deltal_i[i]-DeltalRMS)**2)
		return objective



	def adsdstep(self,forces):
		"""Perform one Adaptive Displacement Steepest Descent step
		Use self.dt as the stepsize
		@param forces: forces array
		"""
		#first calculate the forces directions
		forcesdir=[]
		# since wo store separate geometries per image
		# we must weigh the displacement per image by their
		# total forces. 
		imgforces=[]
		totalforce=0
		for i in forces:
			norm=numpy.dot(i.ravel(),i.ravel())
			totalforce+=norm
			norm=numpy.sqrt(norm)
			imgforces.append(norm)
			forcesdir.append(i.ravel()/norm)
		totalforce=numpy.sqrt(totalforce)
		# now displace the mobile images
		geoshape=numpy.shape(self.geos[0].Geometry)
		for i in range(1,self.numimages()-1):
			self.geos[i].setcoordinates(numpy.reshape((self.geos[i].Geometry.ravel()
				+(self.dt*(imgforces[i]/totalforce))*forcesdir[i]),geoshape))
		# kill possibly stored spline representation
		self.splineRep=None



	def adsdadapt(self,forces):
		"""stepsize adaption procedure for standalone adsd
		@param forces: forces array
		"""
		# can only work, if we have an old force vector to compare with
		stepsizechange=(1.0+(0.618033988/3)) #a lot less than the golden section
		if self.__forces_history_available and self.dtadapted==0:
			fdircheck=0.
			for i in range(1,self.numimages()-1): #ignore fixed images
				fdircheck+=numpy.dot(forces[i].ravel(),self.oldforces[i].ravel())
		else: 
			fdircheck=1.
		if fdircheck<0:
			if self.verbosity>=constants.VBL_NORMAL:
				print("Force reversal detected, reducing stepsize")
			self.dt/=stepsizechange #golden section
			if self.verbosity>=constants.VBL_NORMAL:
				print("New stepsize is {0:12.6f}".format(self.dt))
			self.dtadapted=1
		elif fdircheck > 0.8 and fdircheck < 1.:
			if self.verbosity>=constants.VBL_NORMAL:
				print("Low force direction change, increasing stepsize")
			self.dt *=stepsizechange #goldener schnitt
			if self.verbosity>=constants.VBL_NORMAL:
				print("New stepsize is {0:12.6f}".format(self.dt))
			# paranoia setting: still adapt stepsize in case of overshoot
			self.dtadapted=0
		self.oldforces=forces
		self.__forces_history_available=True



	def adsd(self,forces):
		"""standalone Adaptive Displacement Steepest Descent relaxation
		@param forces: forces array
		"""
		self.adsdstep(forces)
		self.adsdadapt(forces)



	def nebstep(self):
		"""general method to perform one single NEB relaxation step
		using the supplied force calculation and relaxation methods"""
		self.dtadapted=0
		self.nebforces=self.forcefunc()
		self.relaxor(self._fixatoms(self.nebforces))
		# skip overshoot detection when using weighted tangents to avoid
		# spurious geometry resets due to spring force discontinuities
		if self.tangfunc!=self._calcweightedtangents:
			self._checkovershoot(self._fixatoms(self.nebforces))
		self.writecheckpoint(self.checkpointdir)
		self.writexyzpath("lastpath.xyz")
		self.nstep+=1
		# kill possibly stored spline representation
		self.splineRep=None



	def __writePathDebugInfo(self, index):
		"""Write debig information with iteration index "index" to disk
		@param index: iteration number to use as index for debug info
		"""
		global DebugInfoEnergiesFile
		if self._writePathDebugInfoEnergiesFile==None:
			if self.verbosity>=constants.VBL_NORMAL:
				print("writing path debug information to disk")
			DebugInfoEnergiesFile=open("pathDebug.nrg","w")
			self._writePathDebugInfoEnergiesFile=True
			print("#idx     E(max)        Freal(rms)    Freal(max)    Fneb(rms)     Fneb(max)",file=DebugInfoEnergiesFile)
		pathFileNamePrefix="debugpath.{0:04d}".format(index)
		self.writexyzpath(pathFileNamePrefix+".xyz")
		self.writefmgpath(pathFileNamePrefix+".fmg")
		if len(self.nebforces)==0:
			print("{0:4d}  {1: 24.17E}  {2: 24.17E}  {3: 24.17E}  {4: 24.17E}  {5: 24.17E}".format(
				index, max(self.energies), self.rmsforce(force=self.realforces), self.maxforce(force=self.realforces), 0, 0),file=DebugInfoEnergiesFile)
		else:
			print("{0:4d}  {1: 24.17E}  {2: 24.17E}  {3: 24.17E}  {4: 24.17E}  {5: 24.17E}".format(
				index, max(self.energies), self.rmsforce(force=self.realforces), self.maxforce(force=self.realforces), self.rmsforce(force=self.nebforces), self.maxforce(force=self.nebforces)))



	def _sconvcheck(self):
		"""Check Convergence of NEB in standard mode"""
		#We want to use NEB forces on the NEW geometry for convergence checks
		#  but ignore fixed atoms and images
		converged=0
		rms=self.rmsforce(force=self._fixatoms(self.stdnebforces()), images=self.mobrng)
		maxE=max(self.energies[1:self.numimages()-1])
		maxF=self.maxforce(force=self._fixatoms(self.stdnebforces()), images=self.mobrng)
		maxNormal=self.maxforce(force=self._fixatoms(self.realnormalforces()), images=self.mobrng)
		rmsNormal=self.rmsforce(force=self._fixatoms(self.realnormalforces()), images=self.mobrng)
		if self.verbosity>=constants.VBL_NORMAL:
			print("Path RMS NEB force:    {0:12.6f} a.u.\nPath max NEB force:    {1:12.6f} a.u.".format(rms,maxF))
			print("Path RMS normal force: {0:12.6f} a.u.\nPath max normal force: {1:12.6f} a.u.".format(rmsNormal,maxNormal))
			print("Path max Engergy:      {0:12.6f} H".format(maxE))
		if self.energies[0] > self.energies[self.numimages()-1]:
			barrier=abs(maxE-self.energies[0])
		else:
			barrier=abs(maxE-self.energies[self.numimages()-1])
		if self.verbosity>=constants.VBL_NORMAL:
			print("Barrier:               {0:12.6f} H".format(barrier))
		if rms<self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print("path RMS force converged")
			converged=1
		elif  maxF<self.forcetol:
			if self.verbosity>=constants.VBL_QUIET:
				print("path maximum atomic force converged")
			converged=1
		elif self.dt < self.minstepsize:
			if self.verbosity>=constants.VBL_QUIET:
				print("stepsize below {0:6e} -> displacement converged".format(self.minstepsize))
			converged=1
		return converged



	def _ciconvcheck(self):
		"""Check Convergence of NEB in climbing image mode"""
		#We want to use NEB forces on the NEW geometry for convergence checks
		#  but ignore fixed atoms and images
		converged=0
		rms=self.rmsforce(force=self._fixatoms(self.cinebforces()), images=self.mobrng)
		maxF=self.maxforce(force=self._fixatoms(self.cinebforces()), images=self.mobrng)
		maxNormal=self.maxforce(force=self._fixatoms(self.realnormalforces()), images=self.mobrng)
		rmsNormal=self.rmsforce(force=self._fixatoms(self.realnormalforces()), images=self.mobrng)
		CmaxNormal=self.maxforce(force=self._fixatoms(self.realnormalforces()), images=[self.climber])
		CrmsNormal=self.rmsforce(force=self._fixatoms(self.realnormalforces()), images=[self.climber])
		Crms=self.rmsforce(force=self._fixatoms(self.cinebforces()), images=[self.climber])
		CmaxF=self.maxforce(force=self._fixatoms(self.cinebforces()), images=[self.climber])
		if self.verbosity>=constants.VBL_NORMAL:
			print("Path RMS NEB force:     {0:12.6f} a.u.\nPath max NEB force:     {1:12.6f} a.u.".format(rms,maxF))
			print("CI RMS NEB force:       {0:12.6f} a.u.\nCI max NEB force:       {1:12.6f} a.u.".format(Crms,CmaxF))
			print("Path RMS normal force:  {0:12.6f} a.u.\nPath max normal force:  {1:12.6f} a.u.".format(rmsNormal,maxNormal))
			print("CI RMS normal force:    {0:12.6f} a.u.\nCI max normal force:    {1:12.6f} a.u.".format(CrmsNormal,CmaxNormal))
			print("Climbing Image Engergy: {0:12.6f} H".format(self.energies[self.climber]))
		if self.energies[0] > self.energies[self.numimages()-1]:
			barrier=abs(self.energies[self.climber]-self.energies[0])
		else:
			barrier=abs(self.energies[self.climber]-self.energies[self.numimages()-1])
		if self.verbosity>=constants.VBL_NORMAL:	
			print("Barrier:                {0:12.6f} H".format(barrier))
		if Crms<self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print("Climbing image RMS force converged")
			converged=1
		elif  CmaxF<self.forcetol:
			if self.verbosity>=constants.VBL_QUIET:
				print("Climbing image maximum atomic force converged")
			converged=1
		elif self.dt < self.minstepsize:
			if self.verbosity>=constants.VBL_QUIET:
				print("stepsize below {0:8e} -> displacement converged".format(self.minstepsize))
			converged=1
		return converged


	def _checkovershoot(self,forces):
		"""if maximum force or rms force is more than doubled
		by the last step, revert by two geometries and reduce
		stepsize
		@param forces: forces array
		"""
		# only makes sense, if forces and geometries histories are available
		if self.__ovr_history_available and self.__ovr_intermediate_available:
			oldmaxf=self.maxforce(force=self.ovr_intermedforces)
			oldrmsf=self.rmsforce(force=self.ovr_intermedforces)
			####################################################
			#if maximum or RMS doubles between two subsequent steps,
			#revert by two geometries, clear force and geometry histories
			#clear velocities and reduce stepsize
			#######################################################
			if self.maxforce(force=forces) > 2*oldmaxf or self.rmsforce(force=forces) > 2*oldrmsf:
				if self.verbosity>=constants.VBL_NORMAL:
					print("Strong force increase (probably overshoot), reverting geometries and forces")
					#clear histoy
				self.geos=self.ovr_oldgeos
				self.forces=self.ovr_oldforces
				self.__ovr_history_available=False;
				self.__ovr_intermediate_availabe=False;
				self.ovr_intermedforces=None;
				self.ovr_intermedgeos=None;
				self.ovr_oldgeos=None;
				self.ovr_oldforces=None;
				#clear velocities
				for i in range(self.numimages()):
					self.velocities[i]=numpy.zeros((3*self.Atomcount),dtype=float)
				# skip stepzize reduction, if already performed in this step
				if self.dtadapted==0:
					self.dt /= 1.618033988 #goldener schnitt
					self.dtadapted=1
					if self.verbosity>=constants.VBL_NORMAL:
						print("Reducing stepsize due to probable overshoot")
						print("New stepsize = {0:8e}".format(self.dt))
		#store history information
		if self.__ovr_intermediate_available:
			self.ovr_oldgeos = self.ovr_intermedgeos
			self.ovr_oldforces = self.ovr_intermedforces
			self.__ovr_history_available=True
		self.ovr_intermedforces=forces
		self.ovr_intermedgeos=copy.deepcopy(self.geos)
		self.__ovr_intermediate_available=True



	def cubicfit(self):
		"""Fit a cubic polinomial to the path as described in the 
		appendix of J. Chem Phys. 113, 9978 (2000)"""
		#self._calcweightedtangents()
		self._calctangents()
		lengths=[]
		# first calculate path segment lengths, approximated by
		# linear difference between images, projected onto the 
		# vector connecting start- and endstructure
		for i in range(self.numimages()-1):
			xplus=self.geos[i+1].Geometry.ravel()
			x=self.geos[i].Geometry.ravel()
			diff=xplus-x
			R=numpy.dot(diff,diff)
			R=numpy.sqrt(R)
			lengths.append(R)
		# secondly, get the absolute tangent force for each image
		F=[]
		tangforce=self.realtangentforces()
		for i in range (self.numimages()):
			#F.append(sqrt(dot(tangforce[i].ravel(),tangforce[i].ravel())))
			F.append(numpy.dot(tangforce[i].ravel(),self.tangents[i].ravel()))
		# now calculate the Fit constants for each segment
		cubicparms=[]
		for i in range(self.numimages()-1):
			a=(-2*self.energies[i+1]+2*self.energies[i])/lengths[i]**3
			a-=(F[i]+F[i+1])/lengths[i]**2

			b=(3*self.energies[i+1]-3*self.energies[i])/lengths[i]**2
			b+=(2*F[i]+F[i+1])/lengths[i]

			cubicparms.append([a,b,-F[i],self.energies[i],lengths[i]])
		return cubicparms



	def cubicInterpolate(self, steps):
		"""Return a cubic interpolated Path object after J. Chem Phys. 113, 9978 (2000)
		@param steps: number of images in the new path
		"""
		# sanity check: forces must be present!
		if not self.has_realforces():
			raise "cubic interpolation attempted on path without forces."
		#calculate tangents and initialize lengths
		self.tangfunc()
		lengths=[]
		# first calculate path segment lengths, approximated by
		# linear difference between images
		for i in range(self.numimages()-1):
			xplus=self.geos[i+1].Geometry.ravel()
			x=self.geos[i].Geometry.ravel()
			diff=xplus-x
			R=numpy.dot(diff,diff)
			R=numpy.sqrt(R)
			lengths.append(R)
		# calculate total length and distance between new images
		totallength=0.
		for i in lengths:
			totallength+=i
		steplength=totallength/(float(steps)-1.)
		# generate an auxiliary list of pathlengths by segment index to facilitate search for the correct parameter set during interpolation
		pathlengths=[0.]
		lastpathlength=0.
		for i in lengths:
			lastpathlength+=i
			pathlengths.append(lastpathlength)
		# get the flattened tangent force vector for each image
		F=self.realtangentforces()
		for i in range(len(F)):
			F[i]=F[i].ravel()
		# now calculate the Fit coefficient vectors for each segment
		cubicparms=[]
		for i in range(self.numimages()-1):
			a=(-2*self.geos[i+1].Geometry.ravel()+2*self.geos[i].Geometry.ravel())/lengths[i]**3
			a-=(F[i]+F[i+1])/lengths[i]**2
			b=(3*self.geos[i+1].Geometry.ravel()-3*self.geos[i].Geometry.ravel())/lengths[i]**2
			b+=(2*F[i]+F[i+1])/lengths[i]
			cubicparms.append([a,b,-F[i]])
		# instatiate return path object
		returnpath=self.__class__(self.checkpointdir,self.fixedatoms,
		self.dt,"s","s",-1,self.springk,
		"d",self.forcetol,self.rmstol,self.maxit,"s",self.charge,self.verbosity)
		# iterate through new image count and populate new gemetry array
		newgeoarray=[]
		newgeoarray.append(self.geos[0]) #Keep start geomerty untouched!
		segment=0 # keep track of the last parametrization segment
		for i in range(1,steps-1): #account for fixed start- and end geometry
			#find parametrization segment to use
			if segment < len(pathlengths)-1:  # check is unnecessary if we are in the last segment and check would fail!
				while (i*steplength>=pathlengths[segment+1]):
					segment+=1
			# calculate parameter offset in segment
			x=(i*steplength)-pathlengths[segment]
			# calculate new coordinates:
			geoshape=numpy.shape(self.geos[segment].Geometry)
			newcoordinates=self.geos[segment].Geometry.ravel()
			print("{0}".format(x))
			newcoordinates+=x*cubicparms[segment][2]
			newcoordinates+=x*x*cubicparms[segment][1]
			newcoordinates+=x*x*x*cubicparms[segment][0]
			newcoordinates=numpy.reshape(newcoordinates,geoshape)
			# store new coordinates and append to new geometres array
			newgeoarray.append(geometry.Geometry(self.geos[segment].Mode, self.geos[segment].Atomcount, self.geos[segment].AtomTypes, self.geos[segment].Origin, self.geos[segment].Lattice, newcoordinates, self.geos[segment].AtomLayers, self.geos[segment].LayerDict, self.geos[segment].AtomCharges, self.geos[segment].AtomSubTypes, self.geos[segment].LPops))
			newgeoarray[-1].writexyz("dbg-{0:05d}.xyz".format(i))
		# append end geometry
		newgeoarray.append(self.geos[-1])
		returnpath.geos=newgeoarray
		returnpath.Atomcount=newgeoarray[0].Atomcount
		returnpath.amasses=newgeoarray[0].getmasses()
		# finally return interpolated path
		return returnpath



	def cubicenergies(self, cubicparms, steps=100, filename="cubic.nrg"):
		"""Print a list of fitted energies 
		containing steps elements into filename
		@param cubicparms: cubic fit parameters as calculated by cubicfit()
		@param steps: total cound of energy values to calculate (default 100)
		@param filename: output filename, default "cubic.nrg"
		"""
		segsteps=int(steps/self.numimages()-1)+1
		CEfile=open(filename,"w")
		xoffset=0 # this stores the displacement covered between earlier steps
		totallength=0 # total length of the path
		for i in cubicparms:
			totallength+=i[4]
		for i in range(self.numimages()-1):
			for j in range(segsteps):
				x = (float(j)/float(segsteps))*cubicparms[i][4]
				y = cubicparms[i][0]*x*x*x
				y+= cubicparms[i][1]*x*x
				y+= cubicparms[i][2]*x
				y+= cubicparms[i][3]
				print("{0: 24.17E} {1: 24.17E}".format((x+xoffset)/totallength,y),file=CEfile)
			xoffset+=cubicparms[i][4]
		print("\n",file=CEfile)
		x=0
		for i in range(self.numimages()):
			print("{0: 24.17E} {1: 24.17E}".format(x/totallength,self.energies[i]),file=CEfile)
			if i < self.numimages()-1:
				x+=cubicparms[i][4]
		CEfile.close()



	def nebiterate(self,calculator):
		"""Iterate NEB until convergence
		@param calculator: calculator function
		"""
		stopsignal=None
		if self.verbosity>=constants.VBL_NORMAL:
			print("Initializing NEB path search.")
		#open energy trajectory output file
		energyfile=open("energies.dat","a")
		#initialize relaxation
		self.initveloverlet()
		#lastmaxE=0
		#get range of mobile images
		self.mobrng = range(1,self.numimages()-1)
		#Initial calculations
		if self.verbosity>=constants.VBL_NORMAL:
			print("{0:d} Images in path, {1:d} mobile images" .format(self.numimages(),len(self.mobrng)))
			print("Initial Energies and forces calculation")
		self.realforcesfunc(calculator,charge=self.charge)
		rms=self.rmsforce(force=self.realforces)
		maxE=max(self.energies)
		maxF=self.maxforce(force=self.realforces)
		self._calctangents()
		maxNormal=self.maxforce(force=self.realnormalforces())
		rmsNormal=self.rmsforce(force=self.realnormalforces())
		if self.verbosity>=constants.VBL_NORMAL:
			print("Initial path RMS real force  : {0:12.6f}\nInitial path max real force  : {1:12.6f}\nInitial path max Engergy: {2:16.6f}".format(rms,maxF,maxE))
			print("Initial path RMS normal force: {0:12.6f}\nInitial path max normal force: {1:12.6f}\n".format(rmsNormal,maxNormal))
		if self.energies[0] > self.energies[self.numimages()-1]:
			barrier=abs(maxE-self.energies[0])
		else:
			barrier=abs(maxE-self.energies[self.numimages()-1])
		if self.verbosity>=constants.VBL_NORMAL:
			print("Initial barrier:{0:12.6f} H".format(barrier))
		for i in range(self.numimages()):
			print("{0:5d}  {1: 24.17E}".format(i,self.energies[i]),file=energyfile)
		print("\n",file=energyfile)
		energyfile.flush()
		#Do some sanity checking on start- and endpoints and warn the user if necessary
		srms=self.rmsforce([0],self._fixatoms(self.realforces))
		if srms>self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print(" *** WARNING: start image rms force {0:12.6f} a.u. larger than tolerace ***".format(srms))
		smaxF=self.maxforce([0],self._fixatoms(self.realforces))
		if smaxF>self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print(" *** WARNING: start image max force {0:12.6f} a.u. larger than tolerace ***".format(smaxF))
		erms=self.rmsforce([self.numimages()-1],self.realforces)
		if erms>self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print(" *** WARNING:   end image rms force {0:12.6f} a.u. larger than tolerace ***".format(erms))
		emaxF=self.maxforce([self.numimages()-1],self.realforces)
		if emaxF>self.rmstol:
			if self.verbosity>=constants.VBL_QUIET:
				print(" *** WARNING:   end image max force {0:12.6f} a.u. larger than tolerace ***".format(emaxF))
		# Now start cycling
		if self.verbosity>=constants.VBL_NORMAL:
			print("Writing inital checkpoint")
		self.writecheckpoint(self.checkpointdir)
##		self.writefmgpath()
		if self.verbosity>=constants.VBL_NORMAL:
			print("Starting NEB iterations")
		for j in range(self.maxit): #@UnusedVariable
			if self.verbosity>=constants.VBL_NORMAL:
				print("NEB Relaxation step: {0:5d}".format(self.nstep))
			# write debug output, if requested
			if self.writePathDebugInfo:
				self.__writePathDebugInfo(self.nstep)
			self.nebstep()
			self.realforcesfunc(calculator,charge=self.charge)
			for i in range(self.numimages()):
				print("{0:5d}  {1: 24.17E}".format(i,self.energies[i]),file=energyfile)
			print("\n",file=energyfile)
			energyfile.flush()
			if self.convcheckfunc()!=0:
				if self.verbosity>=constants.VBL_QUIET:
					print("Stopping due to convergence.")
				break
			elif os.path.isfile("STOP_comatsci"):
				if self.verbosity>=constants.VBL_QUIET:
					print("Stopping due to stopfile.")
				break
			elif self.stopsignal!=None:
				if self.verbosity>=constants.VBL_QUIET:
					print("Stopping due to Signal: {0}.".format(stopsignal))
				break
			elif self.nstep>=(self.maxit):
				if self.verbosity>=constants.VBL_QUIET:
					print("Maximum number of iterations reached. Stopping.")
				break
			sys.stdout.flush()

		#clean up
		energyfile.close()
##		calc.remove_workdir()

		#write trajectory as xyz for convenience
		self.writexyzpath("path.xyz")

		if self.verbosity>=constants.VBL_QUIET:
			print("finished.")



##	def readcheckpoint(self, checkpointdir):
##		"""read checkpoint from directory and prepare NEBpath to continue path search"""
##		Reactionpath.readcheckpoint(self, checkpointdir)



##	def readfmgpath(self, name):
##		"""read path checkpoint from .fmg file and prepare NEBpath to continue path search"""
##		Reactionpath.readfmgpath(self, name)
