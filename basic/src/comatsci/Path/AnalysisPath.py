## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# AnalysisPath.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################
from __future__ import print_function

from comatsci import constants

from ReactionPath import Reactionpath

import numpy

import math

class AnalysisPath(Reactionpath):
	"""Reaction path with additional analysis functionality, for use in fitrep and analysis tools<br>
	Sets some Reactionpath properties that are meaningless for analysis to hardcoded default values"""



	def __init__(self, icmode,charge=0.0,verbosity=constants.VBL_NORMAL):
		"""Initialize AnalysisPath object
		@param: icmode	calculation mode, c.f. Rectionpath doc
		@param: charge=0.0	system total charge to pass to the claculator
		@param: verbosity=constants.VBL_NORMAL	Verbosity level, c.f. constants doc"""
		# call base class init setting dummy values for useless properties (c.f. Reactionpath doc)
		# keep the base class init quiet, we don't want to print our dummy numbers to the user
		Reactionpath.__init__(self, "", "", icmode, 0.00001, 0.00001, 1, charge,constants.VBL_SILENCE)
		# now set our own verbosity
		self._verbosity=verbosity
		#initialization done



	def forceSquares(self, targetpath,opt=None):
		"""return sum of squares of force difference btw self and targetpath
		@param: targetpath:	second path to compare self against
		@param opt: dummy parameter required by API (default None)		
"""
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating forceSquares")
		Ff=numpy.array(self.realforces).ravel()
		Ft=numpy.array(targetpath.realforces).ravel()
		dF=abs(Ff-Ft)
#		dF=reshape(abs(Ff-Ft),(-1,3))
#		dF=MLab.max(dF,1)
#		dF/=MLab.max(reshape(abs(Ft),(-1,3)),1)
		return numpy.dot(dF,dF)



	def forceRms(self, targetpath,opt=None):
		"""return rms of force deifference btw. self and targetpath
		@param: targetpath:	second path to compare self against
		@param opt: dummy parameter required by API (default None)		
"""
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating forceRms")
		if not targetpath.has_forces():
			raise("AnalysisPath: target path has no forces!")
		if not self.has_forces():
			raise("AnalysisPath: AnalysisPath has no forces!")
		return math.sqrt(self.forceSquares(targetpath)/(self.geos[0].Atomcount*self.numImages))



	def energySquares(self, targetpath):
		"""return squared energy difference btw. self and targetpath
		@param: targetpath:	second path to compare self against"""
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating energySquares")
		Ef=numpy.array(self.energies).ravel()
		Et=numpy.array(targetpath.energies).ravel()
		dE=Ef-Et
		dE/=Et
		return numpy.dot(dE,dE)



	def energyRms(self, targetpath):
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating energyRms")
		if not targetpath.has_energies():
			raise("AnalysisPath: target path has no energies!")
		if not self.has_energies():
			raise("AnalysisPath: AnalysisPath has no energies!")
		"""return rms energy difference btw. self and targetpath
		@param: targetpath:	second path to compare self against"""
		return math.sqrt(self.energySquares(targetpath)/(self.numImages))



	def deltaE(self, educt=0, product=-1,targetpath=None):
		"""return energy difference between educt and product image:
		deltaE=(E_product - E_educt)-reference
		default behavior is to use the first image as educt and the last image as product
		@param educt: educt image index (default 0)
		@param: product=-1 product image index
		@param: targetpath=None	if ! None, compare deltaE to targetpath and subtract targetpath.deltaE(educt, product) (default )		
"""
		if targetpath!=None:
			if not targetpath.has_energies():
				raise "AnalysisPath: target path has no energies!"
			if not self.has_energies():
				raise "AnalysisPath: AnalysisPath has no energies!"
			reference=targetpath.deltaE(educt,product)
		else:
			reference=0.0
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating deltaE")
		return self.energies[product]-self.energies[educt]-reference



	def deltaESquare(self, educt=0, product=-1,targetpath=None):
		"""return squared energy difference between educt and product image:
		deltaESquare=(E_product - E_educt)**2
		default behavior is to use the first image as educt and the last image as product
		@param educt: educt image index (default 0)
		@param: product=-1 product image index
		@param: targetpath=None	if ! None, compare deltaE to targetpath and subtract targetpath.deltaE(educt, product) (default )		
"""
		if targetpath!=None:
			if not targetpath.has_energies():
				raise("AnalysisPath: target path has no energies!")
			if not self.has_energies():
				raise("AnalysisPath: AnalysisPath has no energies!")
			reference=targetpath.deltaE(educt,product)
		else:
			reference=0.0
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating deltaESquare")
		#rather inline than call deltaE for performance reasons
		return (self.energies[product]-self.energies[educt]-reference)**2



	def allDeltaESquares(self,targetpath=None,opt=None):
		"""return sum of all squared energy differences in AnalysisPath (without doube counting):
		allSeltaESquares= sum<sub>i=1..N,j=i+1..N</sub>(E_i - E_j)**2
		@param: targetpath=None	if ! None, compare each deltaE to targetpath and subtract targetpath.deltaE(educt, product) (default )
		@param opt: dummy parameter required by API (default None)		
"""
		squares=0.0
		if self._verbosity>=constants.VBL_DEBUG2:
			print("Analysispath: calculating all squared Energy differences")
		for i in range (self.numImages):
			for j in range(i+1,self.numImages):
				squares+=self.deltaESquare(i,j,targetpath)
		return squares



	def sqreaction(self, targetpath=None, opt=None):
		"""return the square of one reaction energy (or reaction energy difference) calculated from current path
		@param targetpath: target path to calculate reaction energy difference (default None)
		@param opt: options dictionary, must contain: <ul>
		<li>educts: list of images on edict side</li>
		<li>producs: list of images on product side<li>
		</ul>
		Each image can occur multiple times on each side of the equation"""
		if opt==None:
			raise "options required for sqReactions!"
		# if targetpath is specified, calculate reference energy, else leave it at 0
		refEnergy=0.0
		reactionEnergy=0.0
		if targetpath!=None:
			for i in opt["products"]:
				refEnergy+=targetpath.energies[i]
			for i in opt["educts"]:
				refEnergy-=targetpath.energies[i]
		# now calculate the reaction energy, as above
		for i in opt["products"]:
			reactionEnergy+=self.energies[i]
		for i in opt["educts"]:
			reactionEnergy-=self.energies[i]
		# subtract reference, square and return
		return (reactionEnergy-refEnergy)**2
