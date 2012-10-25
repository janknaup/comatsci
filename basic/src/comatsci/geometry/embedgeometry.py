##############################################################################
# EmbedGeometry.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under a free for non-commericial use, do not distribute basis
# see file LICENSE for details.
##############################################################################

import numpy

from . import AnalysisGeometry
from .. import constants #,  utils

#import os
#import sys
#import copy
#import math
#import bisect

class qmmmGeometry(AnalysisGeometry):
	"""Geometry object that supports diverse QM/MM linking schemes"""
	
	# list of features provided by this class
	__geoFeatures__=("embed", )
	
	def  __init__(self, iMode="C", iAtomcount=0, iAtomTypes=None, iOrigin=None, 
		iLattice=None, iGeometry=None,**kwargs):
		"""construct QM/MM Geometry object. See Geometry.__init__ documentation!"""
		AnalysisGeometry.__init__(self, iMode, iAtomcount, iAtomTypes, iOrigin, 
		iLattice, iGeometry,**kwargs)
	
	
	
	def HCSLinkedGeometry(self, othergeo, stepsfunction=None, progressfunction=None, distscale=1.0):
		"""return Homogeneous Charge Substraction linked Geometry, linking with simple hydrogens between QM zone self and modified pointcharges distribution othergeo.<br>
		In the HCS scheme, a fraction dQ of the QM zone total charge is removed from the MM zone atom for each link into the QM zone. dQ=Qtot(QM)/N(link). To conserve overall neutrality, the corresponding partial charge dQ is applied to each link atom
		@param othergeo: point charge distribution to modify and link to
		@param stepsfunction: function to report maximum progress value to (default None)
		@param progressfunction: function to report current progress value to (default None)
		@param distscale: scale the linker-hostatom distances by this value (default 1.0)
		@return: Joint Geometry containing self and link atoms in layer DEFAULT and the modified othergeo in layer PCHR"""
		#first create a deep copy of othergeo to make sure othergeo itself is not modified
		tempothergeo=self.__class__()
		tempothergeo.appendgeometryflat(othergeo,0)
		#now get link list and generate linkatom geometry in the normal way
		#(use the original othergeo, to be sure to correctly take periodicity into account
		linklist = self.getLinkList(othergeo, stepsfunction, progressfunction)
		linksgeo = self.getLinkAtomGeometry(othergeo, linklist, progressfunction, distscale)
		#calculate dQ
		dQ=-self.totalcharge()/float(linksgeo.Atomcount)
		#iterate through link list and apply charge differences to temporary point charges and link atoms geometries
		for (dummy,modindex) in linklist:
			tempothergeo.AtomCharges[modindex]-=dQ
		for i in range(linksgeo.Atomcount):
			linksgeo.AtomCharges[i]=dQ
		#finally, build unified geometry object and return
		returngeo=self.__class__()
		#  copy some attributes from self
		returngeo.Mode=self.Mode
		returngeo.Lattice=self.Lattice
		returngeo.Origin=self.Origin
		#  self into default (QM Zone)
		returngeo.appendgeometryflat(self,0)
		#  add new layer PCHR for modified point charges
		pchrindex=returngeo.addlayer("PCHR")
		returngeo.appendgeometryflat(tempothergeo,pchrindex)
		#  link atoms into default (QM Zone)
		returngeo.appendgeometryflat(linksgeo,0)
		#  finish
		info={
			"chargeTransfer"  : dQ ,
			"linkCount"       : len(linklist),
		}
		return (returngeo, info)



	def BCTCLinkedGeometry(self, othergeo, stepsfunction=None, progressfunction=None, chargeTransfers=None,neutralize=True, distscale=1.0):
		"""return Bond Charge Transfer  linked Geometry, linking with simple hydrogens between QM zone self and modified pointcharges distribution othergeo.<br>
		In the HCS scheme, a fraction dQ of the QM zone total charge is removed from the MM zone atom for each link into the QM zone. dQ=Qtot(QM)/N(link). To conserve overall neutrality, the corresponding partial charge dQ is applied to each link atom
		@param othergeo: point charge distribution to modify and link to
		@param stepsfunction: function to report maximum progress value to (default None)
		@param progressfunction: function to report current progress value to (default None)
		@param ChargeTransfers: Dictionary of element-element charge transfers. Automatically determined by linear-least-squares fit by default (default None)
		@param Neutralize: Boolean specifiying, whether the linked geometry should be neuralized by an additional homogeneous compensating charge transfer per bond (default True)
		@param distscale: scale the linker-hostatom distances by this value (default 1.0)
		@return: tuple of (Joint Geometry containing self and link atoms in layer DEFAULT and the modified othergeo in layer PCHR) and (dictionary containing additional info)"""
		#first create a deep copy of othergeo to make sure othergeo itself is not modified
		tempothergeo=self.__class__()
		tempothergeo.appendgeometryflat(othergeo,0)
		#now get link list and generate linkatom geometry in the normal way
		#(use the original othergeo, to be sure to correctly take periodicity into account
		linklist = self.getLinkList(othergeo, stepsfunction, progressfunction)
		linksgeo = self.getLinkAtomGeometry(othergeo, linklist, progressfunction, distscale)
		#if not given by caller, determine charge transfers automatically
		if chargeTransfers==None:
			chargeTransfers=self.getElementElementChargeTransfers()[0]
		#iterate through link list and apply charge differences to temporary point charges and link atoms geometries
		for i in range(len(linklist)):
			# charge transfer for this particular link, remember that the charge transfer coefficient is stored with index (a,b), a<=b
			if self.AtomTypes[linklist[i][0]]<tempothergeo.AtomTypes[linklist[i][1]]:
				dq=chargeTransfers[(self.AtomTypes[linklist[i][0]],tempothergeo.AtomTypes[linklist[i][1]])] 
			else:
				dq=-chargeTransfers[(tempothergeo.AtomTypes[linklist[i][1]],self.AtomTypes[linklist[i][0]])] 
			#apply charge transfer to MMHA and QML
			tempothergeo.AtomCharges[linklist[i][1]]+=dq
			linksgeo.AtomCharges[i]-=dq
		# calculate residual charge
		Qres=self.totalcharge()+linksgeo.totalcharge()
		# if neutralization is requested, calculate homogeneous charge transfer coefficient and apply to QMHA and QML
		if neutralize:
			#calculate dQ
			dQ=Qres/float(linksgeo.Atomcount)
			for i in range(len(linklist)):
				tempothergeo.AtomCharges[linklist[i][0]]+=dQ
				linksgeo.AtomCharges[i]-=dQ
		else:
			dQ=0.0
		#finally, build unified geometry object and return
		returngeo=self.__class__()
		#  copy some attributes from self
		returngeo.Mode=self.Mode
		returngeo.Lattice=self.Lattice
		returngeo.Origin=self.Origin
		#  self into default (QM Zone)
		returngeo.appendgeometryflat(self,0)
		#  add new layer PCHR for modified point charges
		pchrindex=returngeo.addlayer("PCHR")
		returngeo.appendgeometryflat(tempothergeo,pchrindex)
		#  link atoms into default (QM Zone)
		returngeo.appendgeometryflat(linksgeo,0)
		#  finish
		info={
			"usedChargeTransfers"  : chargeTransfers,
			"BCTCResidualCharge"   : Qres,
			"homoChargeTransfer"   : dQ,
			"linkCount"            : len(linklist),
			"finalCharge"          : returngeo.totalcharge()
			}
		return (returngeo, info)




	def simple_linkatoms(self, othergeo, stepsfunction=None, progressfunction=None, distscale=1.0):
		"""return a geometry containing link hydrogens between self and othergeo on the self side.
		<em>self must not contain any atoms from othergeo</em> If you want linkatoms between two layes
		of geometry do something like geometry.layersubgeometry(layer1).simple_linkatoms(geometr.layersubgeometry(layer2))
		@param center: the geometry to which the returned link atoms should connect self
		@param stepsfunction: callback function to report the total number of progress steps to. For progress display purposes (default None)
		@param progressfunction: callback function to report actual progress to. For progress display purposes. (default None)
		@param distscale: scale the linker-hostatom distances by this value (default 1.0)
		"""
		linklist = self.getLinkList(othergeo, stepsfunction, progressfunction)
		returngeo = self.getLinkAtomGeometry(othergeo, linklist, progressfunction, distscale)
		return returngeo



	def getLinkAtomGeometry(self, othergeo, linklist, progressfunction=None, distscale=1.0):
		"""Build and return a Geometry, containing Linkatoms between self and othergeo on self side along the link specified in linklist
		@param othergeo: geometry object to link to
		@param linklist: list of tuples of the form (selfindex,otherindex), each representing a bond between self's atom[selfindex] and othergeo's atom[otherindex]
		@param progressfunction: function to report progress to, ignore if==None (default None)
		@param distscale: scale the linker-hostatom distances by this value (default 1.0)
		@return: Geometry object containing the link atoms
		"""
		#now construct and populate the link atom geometry
		returngeo=self.__class__()
		linkprogressstep=35.0/float(len(linklist))
		linkprogress=0.0
		for (selfatom,otheratom) in linklist:
			#first the direction from the linked atom to he linkatom
			linkdir=othergeo.Geometry[otheratom]-self.Geometry[selfatom]
			linkdir/=numpy.sqrt(numpy.dot(linkdir,linkdir))
			#calculate the linkatom position at equilibrium bond distance from linked atom
			#Covalent radii table is in Angstrom!
			linklen=(self.CORAD[self.AtomTypes[selfatom]]+self.CORAD[1])*distscale
			linkpos=self.Geometry[selfatom]+(linkdir*(linklen/constants.ANGSTROM))
			#now add the linkatoom
			returngeo.addatom(1,linkpos,subtype="H_l")
			linkprogress+=linkprogressstep
			if progressfunction!=None:
				progressfunction(30+35+int(linkprogress))
		return returngeo



	def SLALinkedGeometry(self, othergeo, stepsfunction=None, progressfunction=None, distscale=1.0):
		"""return linked Geometry, linking with simple hydrogens between QM zone self and modified pointcharges distribution othergeo.<br>
		No manipulation of charges is performed
		@param othergeo: point charge distribution to modify and link to
		@param stepsfunction: function to report maximum progress value to (default None)
		@param progressfunction: function to report current progress value to (default None)
		@param ChargeTransfers: Dictionary of element-element charge transfers. Automatically determined by linear-least-squares fit by default (default None)
		@param Neutralize: Boolean specifiying, whether the linked geometry should be neuralized by an additional homogeneous compensating charge transfer per bond (default True)
		@param distscale: scale the linker-hostatom distances by this value (default 1.0)
		@return: tuple of (Joint Geometry containing self and link atoms in layer DEFAULT and the modified othergeo in layer PCHR) and (dictionary containing additional info)"""
		#build unified geometry object and return
		returngeo=self.__class__()
		#  copy some attributes from self
		returngeo.Mode=self.Mode
		returngeo.Lattice=self.Lattice
		returngeo.Origin=self.Origin
		#  copy self into default (QM Zone)
		returngeo.appendgeometryflat(self,0)
		# generate link Atoms and append to return geometry
		linkAtomGeo=self.simple_linkatoms(othergeo,stepsfunction=stepsfunction, progressfunction=progressfunction,distscale=distscale)
		returngeo.appendgeometryflat(linkAtomGeo,0)
		#  add new layer PCHR for point charges and copy othergeo there
		pchrindex=returngeo.addlayer("PCHR")
		returngeo.appendgeometryflat(othergeo,pchrindex)
		#  finish
		info={
			"QMZcharge"			: self.totalcharge(),
			"PCHRcharge"			: othergeo.totalcharge(),
			"linkCount"            : linkAtomGeo.Atomcount,
			"finalCharge"          : returngeo.totalcharge()
			}
		return (returngeo, info)
