##############################################################################
# erepcalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
from .calculator import Calculator,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED,CALCSTATUS_ERROR,CALCSTATUS_DISABLED #@UnusedImport

from .calcerror import CalcError
from .. import constants
#import comatsci.utils as utils
#import ConfigParser
#import tempfile
#import os
#import sys
#import shutil
import numpy
#import time	
import copy
from comatsci import Spline

class erepcalc(Calculator):
	"""calculates total energy and force contribution from a given
	DFTB repulsive potential"""



	def __init__(self, Erep, verbosity=0):
		"""Erep must be specified, it is the only sensible option for erepcalc
		Erep must be specified in the form of a dictionary mapping interger tuples(i,j)
		to of arrays[2][N] values, where i,j are element numbers with i<=j
		and the arrays contain total energies in Hartree over radii in Bohr.
		@param Erep: dictionary of Repulsive potentials
		@param verbosity: c.f. base class (default 1)		
"""
		if not type(Erep) is dict:
			raise CalcError("Erep is not of type dictionary")
		else:
			self._setEreps(Erep)
		Calculator.__init__(self,verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("initializing E_rep calculator")
		self._status=CALCSTATUS_READY


	# We don't need to implement most of the template functions
	# since we don't call an external caclulator program
	# We must reimplement start(), finreready(), errorreready()
	# and runfg() though

	def start(self,Geometry, steplabel, charge=0, Erep=None):
		"""Calculte energies and forces. Set new E_rep if specified
		@param Geometry: c.f. base class
		@param steplabel: c.f. base class
		@param charge: ignored (default 0)
		@param Erep: if specified, set repulsive potential dictionary (default None)
		prior to calculation"""
		if Erep!=None:
			self._setEreps(Erep)
		tempforces=numpy.zeros((Geometry.Atomcount,3),dtype=float)
		tempenergy=0.0
		#workaround for possible problems with small supercells
		#if any lattice vector is smaller than the largest E_rep cutoff,
		#the Geometry must be expanded in that direction to include enough
		#periodic images
		orig_atc=Geometry.Atomcount
		workgeo=copy.deepcopy(Geometry)
		if Geometry.Mode=="S":
			cutoff=min(self.erep_outercut)
			expand=[0,0,0]
#			lvls=[]
			for i in range(3):
				lvlen=numpy.sqrt(Geometry.Lattice[i][0]**2
					+Geometry.Lattice[i][1]**2+Geometry.Lattice[i][2]**2)
				expand[i]=int(numpy.ceil(cutoff[0]/lvlen))
			workgeo.periodicexpand(expand)
		#END workaround
		dm=workgeo.distancematrix()
		for i in range(workgeo.Atomcount):
			for j in range (i,workgeo.Atomcount):
				dummy1=workgeo.AtomTypes[i]
				dummy2=workgeo.AtomTypes[j]
				if dummy2>dummy1:
					repindex=(dummy1,dummy2)
				else:
					repindex=(dummy2,dummy1)
				dist=dm[i][j]
				if (dist > self.erep_innercut[repindex] and 
				dist < self.erep_outercut[repindex]):
					tempenergy+=self.erep[repindex].splint(dist)
					forcedir=workgeo.Geometry[i]-workgeo.Geometry[j]
					forcedir/=numpy.sqrt(numpy.dot(forcedir,forcedir))
					force= -self.erep[repindex].splder(dist)
					tempforces[i%orig_atc]+=force*forcedir
					tempforces[j%orig_atc]-=force*forcedir
		self.etot=tempenergy
		self.gradients=tempforces
		self._status=CALCSTATUS_FINISHED



	def runfg(self,Geometry, steplabel, charge=0, Erep=None):
		"""directly calculates energy and forces, see start()"""
		self.start(Geometry, steplabel, charge, Erep)



	def _setEreps(self, indict):
		"""set internal Erep spline representations from sample data in indict
		@param indict: E_rep dictionary to set"""
		self.erep={}
		self.erep_innercut={}
		self.erep_outercut={}
		for i in indict.keys():
			self.erep[i]=Spline.spline(indict[i][0],indict[i][1])
			self.erep_innercut[i]=indict[i][0][0]
			self.erep_outercut[i]=indict[i][0][len(indict[i][0])-1]
