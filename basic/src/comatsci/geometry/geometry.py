## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# Geometry.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
from .. import constants, utils
#from numpy.oldnumeric import *
import geoext as gx #@UnresolvedImport
import h5py
import uuid
import numpy
import numpy.linalg as linalg
import os

##import sys
import copy
import math
##import xml.dom.minidom
##import bisect

#complicated import statement to make it work with python 2.4 and 2.5
#  see if python 2.5's elementTree implementation is present
try:
	from xml.etree import ElementTree as ET  # @UnusedImport
#  otherwise try to import locally installed elementtree (for python 2.4 and below)
except:
	from elementtree import ElementTree as ET #@UnresolvedImport @Reimport


##import numpy.oldnumeric.linear_algebra as lina


# some constants for unit conversion
# Bohr Radius to angstrom
Angstrom=constants.ANGSTROM
# Bohr Radius to Bohr Radius (duh)
Bohr=constants.BOHR
Pi=constants.PI


#
# DTD for fmg File format
#
FMG_DTD="""<!ELEMENT mode (#PCDATA)>
<!ELEMENT lattice (latvec_a, latvec_b, latvec_c)>
<!ATTLIST lattice
    orgx CDATA "0.0"
    orgy CDATA "0.0"
    orgz CDATA "0.0"
    lunit (ang|au) "ang"
    >
<!ELEMENT latvec_a (#PCDATA)>
<!ELEMENT latvec_b (#PCDATA)>
<!ELEMENT latvec_c (#PCDATA)>
<!ELEMENT layer (li, lname)>
<!ELEMENT li (#PCDATA)>
<!ELEMENT lname (#PCDATA)>
<!ELEMENT atom (x, y, z, el, st?, chr?, li?, lpop?)>
<!ATTLIST atom
    lunit (ang|au) "ang"
    >
<!ELEMENT x (#PCDATA)>
<!ELEMENT y (#PCDATA)>
<!ELEMENT z (#PCDATA)>
<!ELEMENT el (#PCDATA)>
<!ELEMENT st (#PCDATA)>
<!ELEMENT chr (#PCDATA)>
<!ELEMENT nrg (#PCDATA)>
<!ELEMENT lpop (#PCDATA)>
<!ATTLIST nrg
	eunit (eV|au) "au"
	>
<!ELEMENT velocities (#PCDATA)>
<!ELEMENT forces (#PCDATA)>
<!ELEMENT name (#PCDATA)>
<!ELEMENT inertia (#PCDATA)>
<!ELEMENT rotationals (#PCDATA)>
<!ELEMENT cart_fmax (#PCDATA)>
<!ELEMENT cart_frms(#PCDATA)>
<!ELEMENT zeropoint_corr (#PCDATA)>
<!ELEMENT E_SCF (#PCDATA)>
<!ELEMENT thermo_corr (#PCDATA)>
<!ELEMENT method (#PCDATA)>
<!ELEMENT basis (#PCDATA)>
<!ELEMENT molid (#PCDATA)>
<!ELEMENT specid (#PCDATA)>
<!ELEMENT InChI (#PCDATA)>
<!ELEMENT SMILES (#PCDATA)>
<!ELEMENT SMARTS (#PCDATA)>
<!ELEMENT molfile (#PCDATA)>
<!ELEMENT location (#PCDATA)>
<!ELEMENT calcresults (method, basis?, nrg?, E_SCF?, zeropoint_corr?, cart_frms?, cart_fmax?, rotationals?, inertia?, location?)>
<!ELEMENT dbinfo (name, molid?, specid?, InChI?, SMILES?, SMARTS?, molfile?, location?)>
<!ELEMENT geometry (mode?, lattice?, layer?, atom+)>
<!ELEMENT trjstep (nrg?,velocities?,forces?)>
<!ELEMENT fmg (geometry+,trjstep*,trjinfo?)>
<!ELEMENT stepcount (#PCDATA)>
<!ELEMENT calculator (#PCDATA)>
<!ELEMENT trjinfo (stepcount?,calculator?)>
"""		



class GeometryError(Exception):
	"""Exception Class for Geometry Objects"""
	pass



class GeoLayer:
	"""Geometry Layer descriptor Class, might be extended in future versions!"""

	def __init__(self, Name=""):
		"""initialize Geometry Layer
		@param Name: layer name (default "")
		"""
		self.Name=Name





class Geometry:
	"""Class to store Molecular Geometries and perform various i/o operations
	G{classtree Geometry}
	"""


	#Periodic Table of Elements
	PTE=["X",
	"H" ,                                                                                "He",
	"Li","Be",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
	"Na","Mg",                                                  "Al","Si","P" ,"S" ,"Cl","Ar",
	"K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
	"Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
	"Cs","Ba",
	"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	          "Lu","Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
	"Fr","Ra",
	"Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
	          "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"]



	#Valence Electrons. semicore states are problematic!, this table ends after Ba!
	VALEL=[0,
	1,2,
	1,2,3,4,5,6,7,8,
	1,2,3,4,5,6,7,8,
	1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
	1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
	1,2]



	#Number of valences per atom !!!EXTREMELY CRUDE!!!, only defined up to 20Ca
	VALENCES=[0,
	1,0,
	1,2,3,4,3,2,1,0,
	1,2,3,4,3,2,1,0,
	1,2]


	#lmax for dftb calculations.semicore states are problematic!, this table ends after Ba!
	LMAX=[0,
	1,1,
	2,2,2,2,2,2,2,2,
	3,3,3,3,3,3,3,3,
	4,4,4,4,4,4,4,4,4,4 ,4 ,4 ,4,4,4,4,4,4,
	5,5,5,5,5,5,5,5,5,5 ,5 ,5 ,5,5,5,5,5,5,
	6,6]



	#atomic masses, rounded to two decimal figures. This table ends after Ba!
	AMASS=[9999.9,
	1.01,	4.00,
	6.94,	9.01,	10.81,	12.01,	14.01,	16.00,	19.00,	20.18,
	23.00,	24.31,	26.98,	28.09,	30.97,	32.07,	35.45,	39.95,
	39.10,	40.08,	44.96,	47.87,	50.94,	52.00,	54.94,	55.85,	58.93,	58.69,	63.55,	65.41,	69.72,	72.64,	74.92,	78.96,	79.90,	83.80,
	85.47,	87.62,	88.91,	91.22,	92.91,	95.94,	98.00,	101.97,	102.91,	106.42,	107.87,	112.41,	114.82,	118.71,	121.76,	127.60,	126.90,	131.29,
	132.91,	137.33]


	#covalent radii, given in Angstroms, rounded to two dec. figures. This table ends afer Ba!
	#VDW radius von noble gases, from webelements.com
	#FIXME:	this should be changed to r_covalent for noble gases as well!!!
	CORAD=[100.0,
	0.37,	1.40,
	1.34,	0.90,	0.82,	0.77,	0.75,	0.73,	0.71,	1.54,
	1.54,	1.30,	1.18,	1.11,	1.06,	1.02,	0.99,	1.88,
	1.96,	1.74,	1.44,	1.36,	1.25,	1.27,	1.39,	1.25,	1.26,	1.21,	1.38,	1.31,	1.26,	1.22,	1.19,	1.16,	1.14,	1.10,
	2.11,	1.92,	1.62,	1.48,	1.37,	1.45,	1.56,	1.26,	1.35,	1.31,	1.53,	1.48,	1.44,	1.41,	1.38,	1.35,	1.33,	2.16,
	2.25,	1.98]
	
	#single bond covalent radii. Table ends after Ba!
	#from P. Pyykko and M. Atsumi, Chem. Eur. J. 15 (2009) 186, doi 10.1002/chem.200800987
	SBCR=[100.0,
	0.32,	0.46,
	1.33,	1.02,	0.85,	0.75,	0.71,	0.63,	0.64,	0.67,
	1.55,	1.39,	1.26,	1.16,	1.11,	1.03,	0.99,	0.96,
	1.96,	1.71,	1.48,	1.36,	1.34,	1.22,	1.19,	1.16,	1.11,	1.10,	1.12,	1.18,	1.24,	1.21,	1.21,	1.16,	1.14,	1.17,
	2.10,	1.85,	1.63,	1.54,	1.47,	1.38,	1.28,	1.25,	1.25,	1.20,	1.28,	1.36,	1.42,	1.40,	1.40,	1.36,	1.33,	1.31,
	2.32,	1.96]


	#Reverste PTE, to get Ordinal Number from Atom Symbol
	#(it's so ugly, 'cuz it was generated automatically from PTE...)
	RPTE={'ge': 32, 'gd': 64, 'ga': 31, 'la': 57, 'li': 3, 'tl': 81, 'tm': 69, 'lr': 103,
	'th': 90, 'ti': 22, 'te': 52, 'tb': 65, 'tc': 43, 'ta': 73, 'yb': 70, 'db': 105, 'dy': 66,
	'xe': 54, 'ds': 110, 'h': 1, 'p': 15, 'x': 0, 'zn': 30, 'eu': 63, 'zr': 40, 'er': 68,
	'ru': 44, 're': 75, 'rf': 104, 'rg': 111, 'ra': 88, 'rb': 37, 'rn': 86, 'rh': 45, 'be': 4,
	'ba': 56, 'bh': 107, 'bi': 83, 'bk': 97, 'br': 35, 'c': 6, 'k': 19, 'o': 8, 's': 16,
	'w': 74, 'os': 76, 'co': 27, 'cm': 96, 'cl': 17, 'ca': 20, 'pa': 91, 'cf': 98, 'ce': 58,
	'cd': 48, 'cs': 55, 'cr': 24, 'cu': 29, 'pr': 59, 'pt': 78, 'pu': 94, 'pb': 82, 'lu': 71,
	'pd': 46, 'po': 84, 'pm': 61, 'hs': 108, 'ho': 67, 'hf': 72, 'hg': 80, 'he': 2, 'md': 101,
	'mg': 12, 'b': 5, 'f': 9, 'mo': 42, 'mn': 25, 'n': 7, 'mt': 109, 'v': 23, 'ac': 89,
	'ag': 47, 'ir': 77, 'am': 95, 'al': 13, 'as': 33, 'ar': 18, 'au': 79, 'at': 85, 'in': 49,
	'ni': 28, 'no': 102, 'na': 11, 'nb': 41, 'nd': 60, 'ne': 10, 'es': 99, 'np': 93, 'fr': 87,
	'sc': 21, 'fe': 26, 'fm': 100, 'i': 53, 'sr': 38, 'kr': 36, 'si': 14, 'u': 92, 'sn': 50,
	'sm': 62, 'y': 39, 'sb': 51, 'sg': 106, 'se': 34}

	# List of molecular geometry features provided by this class
	__geoFeatures__=("base", )

	def __init__(self, iMode="C", iAtomcount=0, iAtomTypes=None, iOrigin=None, 
		iLattice=None, iGeometry=None, iAtomLayers=None, iLayerDict=None, 
		iAtomCharges=None,iAtomSubTypes=None, iLPops=None):
		"""Inititalize Geometry object.
		Known modes are
			- B{c} cluster or molecule
			- B{s} supercell
		@param iMode: Geometry mode 
		@param iAtomcount:  atom count (default 0)
		@param iAtomTypes:  atom types list, orogin, lattice vctors and Geometry numpy.array (default None)
		@param iOrigin: 	supercell origin, always ignored but part of .gen file format (default None)
		@param iLattice:  lattice vectors array (default None)
		@param iGeometry: 	atomic positions array (carthesian bohr) (default None)
		@param iAtomLayers: 	list of atom layer index assigments per atom (default None)
		@param iLayerDict: 	dictionary of (atom layer index: atom layer)-s (default None)
		@param iAtomCharges: 	list of single atom charges. Autoinitialized to zero if None. (default None)
		@param iAtomSubTypes:  list of atom subtype strings. Autoinitialized to element symbols if None (default None)
		@param iLPops:  list of atomic l-shell populations, can be list of lists or None (default None)
		"""
		self.Mode=iMode
		self.Atomcount=iAtomcount
		if iAtomTypes==None:
			self.AtomTypes=[]
		else:
			self.AtomTypes=iAtomTypes
		if iOrigin==None:
			self.Origin=numpy.array((0.,0.,0.),dtype=float)
		else:
			self.Origin=iOrigin
		if iLattice==None:
			self.Lattice=numpy.array(([1.,0.,0.],[0.,1.,0.],[0.,0.,1.]),dtype=float)
		else:
			self.Lattice=iLattice
		#default layer is 0
		if iAtomLayers==None:
			self.AtomLayers=[0 for s in range(self.Atomcount)]  # @UnusedVariable
		else:
			self.AtomLayers=iAtomLayers
		#if no Layers are specified, create default layer
		if iLayerDict==None:
			defaultlayer=GeoLayer("default layer")
			self.LayerDict={0: defaultlayer}
		else:
			self.LayerDict=iLayerDict
		# default atom charge is 0
		if iAtomCharges==None:
			self.AtomCharges=[float(0) for s in range(self.Atomcount)]  # @UnusedVariable
		else:
			self.AtomCharges=iAtomCharges
		if iGeometry==None:
			self.Geometry=numpy.zeros((0,3),dtype=float)
		else:
			self.Geometry=iGeometry
		if iAtomSubTypes==None:
			self.AtomSubTypes=[]
			for i in self.AtomTypes:
				self.AtomSubTypes.append(self.PTE[i])
		else:
			self.AtomSubTypes=iAtomSubTypes
		if iLPops==None:
			self.LPops=[]
			for i in range(self.Atomcount):
				self.LPops.append([])
		else:
			self.LPops=iLPops
		self._elementSubGeos={}
		self.uuid=uuid.uuid4()
		self.label=None
		self.method=None
		self.totalenergy=None
		self.ionkineticenergy=None
		self.iontemperature=None
		self.electrontemperature=None
		self.timestep=None
		self.simtime=None
		self.forces=None
		self._reset_derived()
		self._consistency_check()
	
	
	# declare features:
	@classmethod
	def getFeatures(myclass):
		"""return the list of Geometry class features
		@return: list of class features"""
		# initialize features list
		features=[]
		# get inherited features
		for i in myclass.__bases__:
			features+=(i.getFeatures())
		# append own features
		features+=myclass.__geoFeatures__
		# return full list of features
		return features



	def _consistency_check(self):
		"""check data tables for consistency, raise a GeometryError if inconsistencies are encountered
		@return: 1 if all tables are consistent"""
		# check AtomLayers and LayerList for consistency
		for i in self.AtomLayers:
			if not self.LayerDict.has_key(i):
				raise GeometryError('Atom layer list inconsistency')
		# check atom data tables for consistency
		if self.Atomcount != (len(self.Geometry)):
			raise GeometryError('Coordinate count mismatch')
		elif len(self.AtomTypes) != self.Atomcount:
			raise GeometryError('Atom types count mismatch')
		elif len(self.AtomLayers) != self.Atomcount:
			raise GeometryError('Atom layers count mismatch')
		elif len(self.AtomCharges) != self.Atomcount:
			raise GeometryError('Atom charges count mismatch')
		elif len(self.AtomSubTypes) != self.Atomcount:
			raise GeometryError('Atom subtypes count mismatch')
		elif self.LPops!=None:
			if len(self.LPops)!=self.Atomcount:
				raise GeometryError('Atom l-shell populations mismatch')
		elif self.forces!=None:
			if self.forces.shape != self.Geometry.shape:
				raise GeometryError('Force array does not match coordinates')
		return 1



	def _reset_derived(self):
		"""reset derived data arrays"""
		self._blist=None
		self._dmat=None
		self._blmat=None
		self._lal={}
		self._elementSubGeos={}



	def _layeratoms(self, layer):
		"""return list of atoms contained in layer
		@param layer: Layer index from which to generate subgeometry
		"""
		#reuse lists here, function could be called very often
		if not self._lal.has_key(layer):
			tmplist=[]
			for i in range(self.Atomcount):
				if self.AtomLayers[i]==layer:
					tmplist.append(i)
			self._lal[layer]=tmplist
		return self._lal[layer]



	def layersubgeometry(self, layer, mode=None):
		"""return a Geometry object identical to self but containing only
		atoms from layer
		@param layer: layer to create subgeometry from
		@param mode: Mode of the new Geometry, Autoinit to parent Geometry mode if None (default None)
		"""
		atomlist=self._layeratoms(layer)
		if mode==None:
			mode=self.Mode
		lsubgeo=self.__class__(iMode=mode,iLattice=self.Lattice,iOrigin=self.Origin)
		for i in atomlist:
			lsubgeo.addatom(self.AtomTypes[i],self.Geometry[i],None,
				self.AtomCharges[i],self.AtomSubTypes[i])
		return lsubgeo



	def elementsubgeometry(self, element, mode=None, cache=False):
		"""return a Geometry object identical to self but containing only
		atoms of element *element*
		@param element: element to filter for
		@param mode: Mode of the new Geometry, Autoinit to parent Geometry mode if None (default None)
		"""
		if cache and self._elementSubGeos.has_key(element):
			return self._elementSubGeos[element]
		else:
			if mode==None:
				mode=self.Mode
			typesArray=numpy.array(self.AtomTypes)
			subIndices=(typesArray==element).nonzero()[0]
			chargesArray=numpy.array(self.AtomCharges)
			stArray=numpy.array(self.AtomSubTypes)
			esubgeo=self.__class__(iMode=mode,iLattice=self.Lattice,iOrigin=self.Origin,
								iAtomcount=len(subIndices), iAtomTypes=list(typesArray[subIndices]),
								iGeometry=self.Geometry[subIndices],
								iAtomSubTypes=list(stArray[subIndices]))
			esubgeo.parentmap=dict.fromkeys(list(range(esubgeo.Atomcount)),subIndices)
			if cache:
				self._elementSubGeos[element]=esubgeo
		return esubgeo
	
	
	
	def addatom(self, type, position, layer=None, charge=0.0, subtype=None,LPop=None,checkConsistency=True): #FIXME: rename type parameter
		"""append atom of type with position to the geometry
		@param type: atom type (element number)
		@param position: carthesian bohr coordinates
		@param layer: layer index to add atom to (default layer if omitted) (default None)
		@param charge: atom charge (default 0.0)
		@param subtype: Atom Subtype, element symbol if None (default None)
		@param LPop: Atomic l-shell populations list, can aslo be empty or None (default None)
		@param checkConsistency: perform a consistency check of the geometry after appending the atom (default True)
		"""
		if not (layer in self.LayerDict or layer==None):
			raise GeometryError("Trying to add to nonexistent layer")
		self._reset_derived()
		self.AtomTypes.append(type)
		#tempgeo=numpy.zeros((self.Atomcount+1,3),dtype=float)
		tempgeo=numpy.empty((self.Atomcount+1,3),dtype=float)
		if self.Atomcount > 0:
			tempgeo[0:self.Atomcount]=self.Geometry
		tempgeo[self.Atomcount]=numpy.array(position)
		self.Atomcount+=1
		self.Geometry=tempgeo
		if layer==None:
			self.AtomLayers.append(0)
		else:
			self.AtomLayers.append(layer)
		self.AtomCharges.append(float(charge))
		if subtype==None:
			subtype=self.PTE[type]
		self.AtomSubTypes.append(subtype)
		if LPop==None:
			self.LPops.append([])
		else:
			self.LPops.append(LPop)
		if checkConsistency: self._consistency_check()



	def delatom(self, atomno):
		"""delete specified atom from the geometry
		@param atomno:  index of atom to delete (starting from zero!)
		"""
		self._reset_derived()
		del self.AtomTypes[atomno]
		del self.AtomCharges[atomno]
		del self.AtomLayers[atomno]
		del self.AtomSubTypes[atomno]
		del self.LPops[atomno]
		tempgeo=numpy.zeros((self.Atomcount-1,3),dtype=float)
		offset=0
		for i in range(self.Atomcount):
			if i!=atomno:
				tempgeo[i-offset]=self.Geometry[i]
			else:
				offset=1
		self.Geometry=tempgeo
		self.Atomcount-=1
		self._consistency_check()



	def addlayer(self, layername, layerindex=None):
		"""add a new layer to the geometry. return index number of new layer
		@param layername: Name of the layer to add
		@param layerindex: Index of Layer to add. Defaults to  (default None)
		max(LayerDict,Keys())+1, raises GeometryError if index already exists
		"""
		if layerindex==None:
			layerindex=max(self.LayerDict.keys())+1
		elif layerindex in self.LayerDict:
			raise GeometryError("Layer already exists")
		self.LayerDict[layerindex]=GeoLayer(layername)
		return layerindex



	def appendgeometryflat(self, appendgeometry, appendlayer=0):
		"""append geometry to self. Ignore layers in appended geometry
		@param appendgeometry: Geometry object to append
		@param appendlayer: layer index to add atoms to (default 0)
		"""
		if not appendlayer in self.LayerDict.keys():
			raise GeometryError("Trying to add to nonexistent layer")
		for i in range(appendgeometry.Atomcount):
			self.addatom(appendgeometry.AtomTypes[i],
						appendgeometry.Geometry[i],
						appendlayer,
						appendgeometry.AtomCharges[i],
						appendgeometry.AtomSubTypes[i])
	
		
	def readgen(self,filename):
		"""
		read Geometry from file in .gen format
		@type filename: string
		@param filename: name of the file to read
		"""
		infile=utils.compressedopen(filename, "r")
		genString="".join(list(infile))
		infile.close()
		self.parseGenString(genString)
		
	
			
	def parseGenString(self,genstring):
		"""
		Parse geometry representation from string in DFTB .gen format
		@type genstring: string
		@param genstring: string containign DFTB .gen format geometry
		"""
		# split into lines, remove leading and trailing linebreaks
		genLines=genstring.strip("\n").split("\n")
		#parse atom count and geometry mode from fist line
		dummy=genLines[0].split()
		try:
			self.Atomcount=int(dummy[0])
			self.Mode=dummy[1]
			self.Mode=self.Mode.strip()
			self.Mode.capitalize()
		except:
			raise GeometryError("Error parsing first line of .gen file")
		if not self.Mode in ('S','C','F'):
			raise GeometryError("Error in gen file: unknown mode {0:s}".format(self.Mode))
		# consistency-check input file
		if ((self.Mode in ('S','F') and len(genLines)!=self.Atomcount+6) or (self.Mode=='C' and len(genLines)!=self.Atomcount+2)):
			raise GeometryError("Number of lines in .gen file does not match number of atoms")
		# parse atom type symbols from second line
		dummy=genLines[1].split()
		AtomSymbols=[ s.lower() for s in dummy ]
		# iterate through Atom lines
		geo=[]
		self.LPops=[]
		self.AtomTypes=[]
		for i in range(self.Atomcount):
			dummy = genLines[i+2].split()
			try:
				geo.append([ float(s) for s in dummy[2:5] ])
			except:
				raise GeometryError("Error parsing Atom position in line {0:d} of gen file".format(i+1))
			try:
				self.AtomTypes.append(int(self.RPTE[AtomSymbols[int(dummy[1])-1]]))
			except:
				raise GeometryError("Unable to parse atom type specification in line {0:d} of gen file.".format(i+1))
			self.LPops.append([])
		# convert Positions vector to numpy array
		self.Geometry=numpy.array((geo))
		# if direct coordinate mode, convert positions from Angstrom to atomic units
		if self.Mode in ('S','C'):
			self.Geometry/=Angstrom
		# if periodic geometry, parse supercell vectors
		if self.Mode in ('S','F'):
			dummy=genLines[self.Atomcount+2].split()
			try:
				self.Origin=numpy.array( [ float(s)/Angstrom for s in dummy[:3] ] )
			except:
				raise GeometryError("Error parsing periodic cell origin in gen file")
			try:
				for i in range(3):
					dummy=genLines[self.Atomcount+3+i].split()
					self.Lattice[i]=numpy.array( [ float(s)/Angstrom for s in dummy[:3] ] )
			except:
				raise GeometryError("Error parsing supercell vectors in gen file")
		# convert from fractional to direct coordinates
		if self.Mode=="F":
			self.Geometry = numpy.dot(self.Geometry,self.Lattice)
			self.Mode="S"
		# add dummy supercell spec to cluster Geometry
		elif self.Mode=="C":
			self.Origin=numpy.array((0.,0.,0.),dtype=float)
			self.Lattice=numpy.array(([1.,0.,0.],[0.,1.,0.],[0.,0.,1.]),dtype=float)
		# supply default datastructures which cannot be read from .gen file
		defaultlayer=GeoLayer("default layer")
		self.LayerDict={0: defaultlayer}
		self.AtomCharges=[float(0) for s in range(self.Atomcount)]
		self.AtomLayers=[0 for s in range(self.Atomcount)]
		self.AtomSubTypes=[self.PTE[self.AtomTypes[s]] for s in range(self.Atomcount)]
		self._consistency_check()



	def getAimsString(self):
		"""
		@rtype: string
		@return: geometry as FHI AIMS geometry.in string 
		""" 
		# initialize list of return lines
		rlines=["# geometry file created by comatsci",]
		for i in range(self.Atomcount):
			# convert coordinates to Angstrom
			coordinates=numpy.array(self.Geometry[i])
			coordinates*=constants.ANGSTROM
			# if AtomSubType is defined, use that as label, otherwise use Element Symbol
			try:
				typelabel=self.AtomSubTypes[i]
			except IndexError:
				typelabel=self.PTE[self.AtomTypes[i]]
			# assemble atom line
			rlines.append("atom   {0[0]: 24.17E}   {0[1]: 24.17E}   {0[2]: 24.17E}   {1:s}".format(coordinates,typelabel))
		# write lattice vectors if supercell geometry
		if self.Mode=="S":
			for i in range(3):
				#convert lattice vector to Angstrom
				lv=numpy.array(self.Lattice[i])
				lv*=constants.ANGSTROM
				# assemble lattice line
				rlines.append("lattice_vector   {0: 24.17E}   {1: 24.17E}   {2: 24.17E}   ".format(lv[0],lv[1],lv[2]))
		# finished, return
		return "\n".join(rlines)
	
	
	aimsString=property(fget=getAimsString,fset=None,fdel=None,doc="Geometry string representation in FHI AIMS format")
	
	
	
	def writeAIMS(self,filename="geometry.in"):
		"""
		Write geometry as FHI AIMS inout file
		"""
		outfile=open(filename, "w")
		print(self.aimsString,file=outfile)
		outfile.close()


	knownCDHFields=("uuid","method","totalenergy","ionkineticenergy","iontemperature",
				"electrontemperature","timestep","simtime","coordinates","forces",
				"elements","types","charges","lattice","residues")
	def writeCDHFrameGroup(self,h5file, groupname="frame0000000000",overwrite=False, 
						labelstring="comatsci geometry",exclude=None):
		"""
		write HDF5 representation of Geometry into HDF5 file followin CDH specification
		DOES NOT CLOSE the HDF5 file under any circumstances to allow writing multi-frame geometries into one file
		@type h5file: h5py file object
		@param h5file: HDF5 file to write Geometry into. If no file is passed, a new file named according to parameter filename will be created
		@type groupname: string
		@param groupname: name of the group to write data into.   
		@type overwrite: boolean
		@param overwrite: if true, overwrite preexisting geometry group
		@type labelstring: string
		@param labelstring: string to be attached to the HDF frame group as the label property
		@type exclude: sequence
		@param exclude: list of attribute and set names to exclude 
		"""
		# check overwrite 
		# *** deactivated to reduce keys() calls 
# 		if groupname in h5file.keys() and overwrite==False:
# 			raise ValueError("specified geometry group already present in HDF5 file")
		# ...
		if exclude==None: exclude=[]
		# create HDF5 datatypes
		ResidueDict={}
		for i in self.LayerDict.keys():
			ResidueDict[self.LayerDict[i].Name]=i
		# hack for h5py API changes
		# *** depending on installed h5py version, one of this will cause errors, hence UndefinedVariable comments ***
		if h5py.version.version_tuple[0]==1 and h5py.version.version_tuple[1]<3:
			vlstring=h5py.new_vlen(str)#@UndefinedVariable @UnusedVariable
			layerenum=h5py.new_enum("i",ResidueDict)#@UndefinedVariable
		else: 
			vlstring=h5py.special_dtype(vlen=str) #@UndefinedVariable @UnusedVariable
			layerenum=h5py.special_dtype(enum=("i",ResidueDict)) #@UndefinedVariable
		# get frame group object, create if necessary. We already checked collision with existing froup if overwrite==False
		# overwriting datasets will still fail!
		framegroup=h5file.require_group(groupname)
		# create datasets, initialize with data where straightforward
		if labelstring!=None and self.label!=None:
			if self.label!=None:
				framegroup.attrs["label"]=self.label
			else:
				framegroup.attrs["label"]=labelstring
		framegroup.attrs["uuid"]=str(self.uuid)
		if self.method != None and (not "method" in exclude): 
			framegroup.attrs["method"]=self.method
		if self.totalenergy != None and (not "totalenergy" in exclude): 
			framegroup.attrs["totalenergy"]=self.totalenergy
		if self.ionkineticenergy != None and (not "ionkineticenergy" in exclude): 
			framegroup.attrs["ionkineticenergy"]=self.ionkineticenergy
		if self.iontemperature != None and (not "iontemperature" in exclude): 
			framegroup.attrs["iontemperature"]=self.iontemperature
		if self.electrontemperature != None and (not "electrontemperature" in exclude): 
			framegroup.attrs["electrontemperature"]=self.electrontemperature
		if self.timestep != None and (not "timestep" in exclude): 
			framegroup.attrs["timestep"]=self.timestep
		if self.simtime != None and (not "simtime" in exclude): 
			framegroup.attrs["simtime"]=self.simtime
		if (not "coordinates" in exclude): 
			geoset=framegroup.create_dataset("coordinates",data=numpy.array(self.Geometry,'=f8'))  # @UnusedVariable
		if self.forces != None and (not "forces" in exclude):
			forceset=imagegroup.create_dataset("forces",data=num.array(self.forces,"=f8")) #@UnusedVariable
		if self.AtomTypes != None and (not "elements" in exclude):
			elementset=framegroup.create_dataset("elements",data=numpy.array(self.AtomTypes,'=u1')) #@UnusedVariable
		if self.AtomSubTypes != None and (not "types" in exclude):
			typeset=framegroup.create_dataset("types",data=numpy.asarray(self.AtomSubTypes,dtype="string_"))  # @UnusedVariable
		# charges are tricky, but we should just leave them out if all charges in Geometry instance are zero anyway
		if numpy.shape(numpy.nonzero(self.AtomCharges))[1]!=0 and (not "charges" in exclude):
			chargeset=framegroup.create_dataset("charges",data=numpy.array(self.AtomCharges,'=f8')) #@UnusedVariable
		if (not "lattice" in exclude):
			if self.Mode=="S":
				latticeset=framegroup.create_dataset("lattice",data=numpy.array(self.Lattice,'=f8'))  # @UnusedVariable
		# layer indices might not be contigous, depending on geometry history
		if (not "residues" in exclude):
			residueset=framegroup.create_dataset("residues",(self.Atomcount,),dtype=layerenum) #@UnusedVariable
			# populate so far uninitialized data sets
			residueset=numpy.array(self.AtomLayers) #@UnusedVariable
		# to be able to add data into the framegroup created, return references to HDF5 file and data group
		return (h5file,framegroup)
	
	
	
	def writeCDH(self,filename):
		"""
		write Geometry as Chemical Data Hierachy file
		@type filename: string
		@param filename: name of the cdh file to write to
		"""
		h5file=h5py.File(filename,"w")
		self.writeCDHFrameGroup(h5file=h5file,overwrite=True)
		h5file.close()



	def parseH5Framegroup(self, framegroup, globalsGroup=None):
		"""
		set own properties from HDF5 framegroup
		@type framegroup: HDF5 data group following CDH convention
		@param framegroup: geometry representation to parse geometry from
		@type globalsGroup: HDF data group following CDH convention
		@param globalsGroup: if provided, data not present in the frame group 
		will be read from the globas group.Intended to reduce storage space requirements of CDH files. 
		@return: reference to self
		"""
		# if globasGroup is not specified, create an empty dictionary as dummy
		if globalsGroup==None: globalsGroup=framegroup
		# CDH compliance check
		framesets=framegroup.keys()
		globalsets=globalsGroup.keys()
		if not (("coordinates" in framesets or "coordinates" in globalsets) and 
			("elements" in framesets) or ("elements" in globalsets)):
			raise ValueError("non CDH compliant data group passed for parsing.")
		# clear self
		self._reset_derived()
		#********************************************************************
		# Frame attributes first
		#********************************************************************
		# uuid
		if "uuid" in framegroup.attrs.keys():
			self.uuid=uuid.UUID(framegroup.attrs["uuid"])
		# label
		if "label" in framegroup.attrs.keys():
			self.label=str(framegroup.attrs["label"])
		# method
		if "method" in framegroup.attrs.keys():
			self.method=str(framegroup.attrs["method"])
		elif "method" in globalsGroup.attrs.keys():
			self.method=str(globalsGroup.attrs["method"])
		# totalenergy
		if "totalenergy" in framegroup.attrs.keys():
			self.totalenergy=float(framegroup.attrs["totalenergy"])
		elif "totalenergy" in globalsGroup.attrs.keys():
			self.totalenergy=float(globalsGroup.attrs["totalenergy"])
		# ionkineticenergy
		if "ionkineticenergy" in framegroup.attrs.keys():
			self.ionkineticenergy=float(framegroup.attrs["ionkineticenergy"])
		elif "ionkineticenergy" in globalsGroup.attrs.keys():
			self.ionkineticenergy=float(globalsGroup.attrs["ionkineticenergy"])
		# iontemperature
		if "iontemperature" in framegroup.attrs.keys():
			self.iontemperature=float(framegroup.attrs["iontemperature"])
		elif "iontemperature" in globalsGroup.attrs.keys():
			self.iontemperature=float(globalsGroup.attrs["iontemperature"])
		# iontargettemperature
		if "iontargettemperature" in framegroup.attrs.keys():
			self.iontargettemperature=float(framegroup.attrs["iontargettemperature"])
		elif "iontargettemperature" in globalsGroup.attrs.keys():
			self.iontargettemperature=float(globalsGroup.attrs["iontargettemperature"])
		# electrontemperature
		if "electrontemperature" in framegroup.attrs.keys():
			self.electrontemperature=float(framegroup.attrs["electrontemperature"])
		elif "electrontemperature" in globalsGroup.attrs.keys():
			self.electrontemperature=float(globalsGroup.attrs["electrontemperature"])
		# timestep
		if "timestep" in framegroup.attrs.keys():
			self.timestep=int(framegroup.attrs["timestep"])
		elif "timestep" in globalsGroup.attrs.keys():
			self.timestep=int(globalsGroup.attrs["timestep"])
		# simtime
		if "simtime" in framegroup.attrs.keys():
			self.simtime=float(frameGroup.attrs["simtime"])
		elif "simtime" in globalsGroup.attrs.keys():
			self.simtime=float(globalsGroup.attrs["simtime"])
		#********************************************************************
		# Now data sets
		#********************************************************************
		# first check if geometry is cluster or periodic and set lattice if necessary:
		if "lattice" in framesets:
			self.Mode="S"
			self.Lattice=framegroup["lattice"].value
		elif "lattice" in globalsets:
			self.Mode="S"
			self.Lattice=globalsGroup["lattice"].value
		else:
			self.Mode="C"
		# set coordinates
		if "coordinates" in framesets:
			coordgroup=framegroup
		else:
			coordgroup=globalsGroup
		self.Geometry=coordgroup["coordinates"].value
		# set elements
		if "elements" in framesets:
			self.AtomTypes=list(framegroup["elements"].value)
		else:
			self.AtomTypes=list(globalsGroup["elements"].value)
		self.Atomcount=len(self.AtomTypes)
		# set atom subtypes
		if "types" in framesets:
			self.AtomSubTypes=list(framegroup["types"].value)
		elif "types" in globalsets:
			self.AtomSubTypes=list(globalsGroup["types"].value)
		else:
			self.AtomSubTypes=[ self.PTE[self.AtomTypes[ii]] for ii in range(self.Atomcount) ] 
		# set atom layers
		residueGroup=None
		if "residues" in framesets:
			residueGroup=framegroup
		elif "residues" in globalsets:
			residueGroup=globalsGroup
		if residueGroup!=None:
			# hack around api changes in h4py
			if h5py.version.version_tuple[0]==1 and h5py.version.version_tuple[1]<3:
				residuesDict=h5py.get_enum(residueGroup["residues"].dtype) #@UndefinedVariable
			else: 
				residuesDict=h5py.check_dtype(enum=residueGroup["residues"].dtype) #@UndefinedVariable
			self.LayerDict={}
			print(residuesDict)
			for ii in residuesDict.keys():
				if self.layerbyname(residuesDict[ii])==None:
					self.addlayer(ii,int(residuesDict[ii]))
			self.AtomLayers=residueGroup["residues"].value
			print(self.AtomLayers)
		else:
			if not self.layerbyname("default layer")==0:
				self.addlayer("default layer", 0)
			self.AtomLayers=list(numpy.zeros((self.Atomcount,)))
		# set atom charges
		if "charges" in framesets:
			self.AtomCharges=list(framegroup["charges"].value)
		elif "charges" in globalsets:
			self.AtomCharges=list(globalsGroup["charges"].value)
		else:
			self.AtomCharges=list(numpy.zeros((self.Atomcount,),dtype=float))
		# set forces
		if "forces" in framesets:
			forcesgroup=framegroup
		else:
			forcesgroup=globalsGroup
		if "forces" in globalsets:
			self.forces=forcesgroup["forces"].value
		else:
			self.forces=None
		# dummy data
		self.LPops=list(numpy.zeros((self.Atomcount,3),dtype=int))
		# Finally, check consistency and return
		self._consistency_check()
		return self



	def readCDHFile(self,filename,groupname=None):
		"""
		read Geometry from Chemical Data Hierachy file
		@type filename: string
		@param filename: name of the cdh file read from
		@type groupname: string
		@param groupname: name of the HDF5 data group to read from, frame with alphabetically lowest label 
		"""
		h5file=h5py.File(filename,"r")
		# fine the aphabatically first item named frame*
		if groupname==None:
			filekeys=h5file.keys()
			filekeys.sort()
			groupname=filekeys.pop(0)
			while (groupname[0:5].lower()!="frame"):
				groupname=filekeys.pop(0)
		self.parseH5Framegroup(h5file[groupname],h5file.get("globals",None))
		h5file.close()



	def parseAimsString(self,instring):
		""" get geometry information from FHI aims string
		@type instring: string
		@param instring: fhi aims formated geometry specification  
		"""
		#initialize temporary lattice, coordinates and atom types
		templattice=[]
		temptypes=[]
		tempcoords=[]
		tcfractional=[]
		#split string into lines
		lines=instring.split("\n")
		#iterate through lines
		for i in range(len(lines)):
			#tokenize lines, split off comments at EOL if any
			tokens=lines[i].strip().split("#")[0].split()
			#ignore empty and comment lines
			if (len(tokens)!=0) and (len(tokens[0])!=0) and (tokens[0][0]!="#"):
				# we only read lines of type atom, atom_frac and lattice_vector
				if tokens[0].lower()=="lattice_vector":
					if len(tokens)!=4:
						raise ValueError("Error in line {0:d} of FHI aims string: lattice vector needs exactly 3 components".format(i+1))
					else:
						templattice.append([tokens[1],tokens[2],tokens[3]])
				# we support carthesian and fractional coordinates
				elif tokens[0].lower() in ["atom","atom_frac"]:
					if len(tokens)!=5:
						raise ValueError("Error in line {0:d} of FHI aims string: atom needs exactly 3 coordinates and atom type".format(i+1))
					else:
						tempcoords.append(numpy.array((tokens[1],tokens[2],tokens[3])))
						if tokens[0]=="atom":
							tcfractional.append(False)
						elif tokens[0]=="atom_frac":
							tcfractional.append(True)
						temptypes.append(self.RPTE[tokens[4].lower()])
				else:
					raise ValueError("Unrecognized keyword {0:s} in line {1:d} of FHI aims string".format(tokens[0],i+1))
		# Handle lattice vectors. We only support non-periodic and 3d periodic geometries
		if len(templattice)<3:
			self.Mode="C"
			if len(templattice)>0:
				print("Warning, less than 3D periodicity ignored!")
		elif len(templattice)==3:
			self.Mode="S"
			self.Lattice=numpy.array(templattice,dtype=float)/Angstrom
		else:
			raise ValueError("More than three lattice vectors in geometry.")
		# iterate through atoms
		for i in range(len(tempcoords)):
			position=numpy.array(tempcoords[i],dtype=float)
			# if coordinates for current atom are fractional, convert to carthesian
			if tcfractional[i]:
				if self.Mode!="S":
					raise ValueError("Fractional coordinates specified for atom in non 3D periodic geometry.")
				else:
					position=numpy.dot(position,self.Lattice)
			else:
				position/=Angstrom
			self.addatom(temptypes[i], position)
		# finished, clean up derived properties
		self._reset_derived()



	def readAimsFile(self,filename="geometry.in"):
		""" read FHI aims geometry from file
		@type filename: string
		@param filename: name of the file to read. Default is FHI default geometry input file name "geometry.in"  
		"""
		infile=utils.compressedopen(filename,autodetect=False)
		inlines=list(infile)
		infile.close()
		instring="\n".join(inlines)
		self.parseAimsString(instring)
		# done.


	def parseXyzString(self, xyzstring):
		"""parse geometry string in cyz format
		@type xyzstring: string
		@param xyzstring: string containing Atomic coordinated in xmol .xyz format
		"""
		# prepare temporary variables
		tempgeo=[]
		tempAtomTypes=[]
		tempLPops=[]
		tempAtomCharges=[]
		# break string into a list of lines
		xyzlines=xyzstring.split('\n')
		# discard trailing empty line if present
		if xyzlines[-1].strip()=="":
			del xyzlines[-1]
		# parse atom count
		line=xyzlines.pop(0)
		try:
			tempAtomCount=(int(line))
		except:
			print("atom count '{0:s}' in xyz string could not be parsed. Abort".format(line))
			raise
		# discard comment line, then check if number of atom lines
		# matches number of atom
		del xyzlines[0]
		if len(xyzlines)!=tempAtomCount:
			raise(ValueError,"Number of atom lines does not match specified atom count in xyz string")
		# iterate trhough atom lines
		for i in range(tempAtomCount):
			line=xyzlines[i].split()
			atsym=line[0].strip().lower()
			try:
				tempAtomTypes.append(self.RPTE[atsym])
			except:
				print("Atom symbol '{0:s}' in line {0:d} of xyz string could not be parsed. Abort.".format(atsym,i+3))
				raise
			tempLPops.append([])
			try:
				tempgeo.append([ float(s)/Angstrom for s in line[1:4] ])
			except:
				print("Atomic coordinates in line {0:d} of xyz string could not be parsed. Abort.".format(i+3))
				raise
			# DFTB specialty: try parsing coulumn 6 as atomic charge, ignore if unsuccessful
			if len(line)>=5:
				try:
					tempAtomCharges.append(-(float(line[4])-self.VALEL[tempAtomTypes[-1]]))
				except:
					tempAtomCharges.append(0.0)
			else:
				tempAtomCharges.append(0.0)
		# clear self before applying parsed geometry
		self._reset_derived()
		# apply parsed geometry and set implicit defaults
		self.Atomcount=tempAtomCount
		self.AtomTypes=tempAtomTypes
		self.LPops=tempLPops
		self.Mode='C'
		self.Origin=numpy.array((0.,0.,0.),dtype=float)
		self.Lattice=numpy.array(([1.,0.,0.],[0.,1.,0.],[0.,0.,1.]),dtype=float)
		self.Geometry=numpy.array((tempgeo))
		defaultlayer=GeoLayer("default layer")
		self.LayerDict={0: defaultlayer}
		self.AtomCharges=tempAtomCharges
		self.AtomLayers=[0 for s in range(self.Atomcount)]
		self.AtomSubTypes=[self.PTE[self.AtomTypes[s]] for s in range(self.Atomcount)]
		# finally, check self for sanity.
		self._consistency_check()
		# finished.
		



	def readxyz(self, filename):
		"""read geometry from xyz file
		@param filename: input file name
		"""
		self._reset_derived()
		infile=utils.compressedopen(filename,"r")
		instring="".join(list(infile))
		infile.close
		try:
			self.parseXyzString(instring)
		except:
			print("Parsing of .xyz file '{0:s}' failed. Abort.".format(filename))
			raise



	def handlegeometry_dom(self, ingeo):
		"""handle geometry element from fmg file in (mini)dom representation
		@param ingeo: (mini)dom geometry Element
		"""
		#constants etc
		lunits={"ang":Angstrom, "au":Bohr}
		#set some default values
		self.Atomcount=0
		self.Origin=numpy.zeros((3),dtype=float)
		self.Lattice=numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],dtype=float)
		#first get the mode
		mode=ingeo.getElementsByTagName("mode")
		if len(mode)==1:
			dummy=str(mode[0].childNodes[0].data.strip().upper())
			if dummy in ["C","S"]:
				self.Mode=dummy
			else:
				raise GeometryError("FMG: Invalid geometry mode")
		else:
			self.Mode="C"
		# now get the lattice information, if mode is "S"
		if self.Mode=="S":
			lattice=ingeo.getElementsByTagName("lattice")[0]
			if lattice.hasAttribute("lunit"):
				lunit=lunits[lattice.getAttribute("lunit").strip()]
			else:
				#default length unit is Angstrom
				lunit=Angstrom
			if lattice.hasAttribute("orgx"):
				self.Origin[0]=float(lattice.getAttribute("orgx"))/lunit
			if lattice.hasAttribute("orgy"):
				self.Origin[1]=float(lattice.getAttribute("orgy"))/lunit
			if lattice.hasAttribute("orgz"):
				self.Origin[2]=float(lattice.getAttribute("orgz"))/lunit
			latv_a=lattice.getElementsByTagName("latvec_a")
			if len(latv_a)==1:
				dummy=latv_a[0].childNodes[0].data.split()
				for i in range(3):
						self.Lattice[0][i]=float(dummy[i])/lunit
			latv_b=lattice.getElementsByTagName("latvec_b")
			if len(latv_b)==1:
				dummy=latv_b[0].childNodes[0].data.split()
				for i in range(3):
						self.Lattice[1][i]=float(dummy[i])/lunit	
			latv_c=lattice.getElementsByTagName("latvec_c")
			if len(latv_c)==1:
				dummy=latv_c[0].childNodes[0].data.split()
				for i in range(3):
						self.Lattice[2][i]=float(dummy[i])/lunit
		#now get the Layers dictionary
		layers=ingeo.getElementsByTagName("layer")
		if len(layers)>0:
			for i in layers:
				layerid=int(i.getElementsByTagName("li")[0].childNodes[0].data)
				name=i.getElementsByTagName("lname")[0].childNodes[0].data.strip()
				self.LayerDict[layerid]=GeoLayer(name)
		#We need al least the default Layer
		else:
			self.LayerDict={0:GeoLayer("default Layer")}
		#Finally handle all atom records
		atoms=ingeo.getElementsByTagName("atom")
		for i in atoms:
			if i.hasAttribute("lunit"):
				lunit=lunits[i.getAttribute("lunit").strip()]
			#required elements
			pos=[0.,0.,0.]
			pos[0]=float(i.getElementsByTagName("x")[0].childNodes[0].data.strip())/lunit
			pos[1]=float(i.getElementsByTagName("y")[0].childNodes[0].data.strip())/lunit
			pos[2]=float(i.getElementsByTagName("z")[0].childNodes[0].data.strip())/lunit
			el=int(i.getElementsByTagName("el")[0].childNodes[0].data.strip())
			#optional elements
			li=i.getElementsByTagName("li")
			if len(li)==1:
				layer=int(li[0].childNodes[0].data)
			else:
				layer=0
			st=i.getElementsByTagName("st")
			if len(st)==1:
				subtype=st[0].childNodes[0].data.strip()
			else:
				subtype=self.PTE[el]
			chrg=i.getElementsByTagName("chr")
			if len(chrg)==1:
				charge=float(chrg[0].childNodes[0].data)
			else:
				charge=0.0
			lpop=i.getElementsByTagName("lpop")
			if len(lpop)==1:
				dummy=lpop[0].childNodes[0].data.strip().split()
				lpopulation=[ float(s) for s in dummy ]
			else:
				lpopulation=None
			self.addatom(el,pos,layer,charge,subtype,lpopulation)
		self._consistency_check()
		self._reset_derived()



##	def readfmg(self, filename):
##		"""read geometry from fmg file, using minidom
##		@param filename: input file name
##		"""
##		dom = xml.dom.minidom.parse(filename)
##		ingeo=dom.getElementsByTagName("geometry")[0]
##		self.handlegeometry_dom(ingeo)
##		# cleanup to save memory
##		ingeo.unlink()



	def readfmg(self, filename):
		"""read geometry from fmg file, using elementtree
		@param filename: input file name
		"""
		# open the fmg file by outside of ET to allow transparent decompression
		fmgfile=utils.compressedopen(filename)
		# parse xml file into element tree
		tree = ET.parse(fmgfile)
		# close the file object asap
		fmgfile.close()
		# handle the element tree
		self.handleGeoElementTree(tree.getroot().find("geometry"))
		# done



	def handleGeoElementTree(self, tree):
		"""handle xml element tree representation of geometry
		@param tree: element tree to parse
		"""
		# some constants
		lunits={"ang":Angstrom, "au":Bohr}
		allowedmodes=("s","c")
		# get all geometry elements from the tree
		geoIterator=tree.getiterator("geometry")
		# bomb out, if we have more than one geometry in the tree passes here!
		if len(geoIterator)>1:
			raise GeometryError("More than one Geometry in .fmg file")
		# get mode
		mode=tree.findtext("mode")
		#   default mode is "C"
		if mode==None:
			self.Mode="C"
		#   allowed modes are C_luster and S_upercell
		elif mode.lower() in allowedmodes:
			self.Mode=mode.upper()
		else:
			raise GeometryError("invalid geometry mode in .fmg file")
		# get lattice if mode is S
		if self.Mode=="S":
			lattice=tree.find("lattice")
			# lattice must be given if mode is S
			if lattice==None:
				raise GeometryError("No lattice vectors given for periodic geometry!")
			# get length unit, default is ang_strom
			lunit=lunits[lattice.get("lunit","ang").lower().strip()]
			# get origin, default is (0.0, 0.0, 0.0)
			orgx=float(lattice.get("orgx","0.0"))/lunit
			orgy=float(lattice.get("orgy","0.0"))/lunit
			orgz=float(lattice.get("orgz","0.0"))/lunit
			self.Origin=numpy.array((orgx,orgy,orgz))
			# get lattice vectors
			lva=[ float(s)/lunit for s in lattice.findtext("latvec_a").split() ]
			lvb=[ float(s)/lunit for s in lattice.findtext("latvec_b").split() ]
			lvc=[ float(s)/lunit for s in lattice.findtext("latvec_c").split() ]
			self.Lattice=numpy.array([lva,lvb,lvc])
			# done parsing lattice
		# get layer specs, create default layer if none is specified
		layers=tree.getiterator("layer")
		if len(layers)>0:
			self.LayerDict={}
			for i in layers:
				layerid=int(i.findtext("li"))
				name=i.findtext("lname").strip()
				self.LayerDict[layerid]=GeoLayer(name)
		#We need at least the default Layer
		else:
			self.LayerDict={0:GeoLayer("default Layer")}
		# Now iterate through atoms
		atoms=tree.getiterator("atom")
		for i in atoms:
			# start with length unit
			lunit=lunits[i.get("lunit","ang").lower().strip()]
			# get the coordinates
			x=float(i.findtext("x").strip())/lunit
			y=float(i.findtext("y").strip())/lunit
			z=float(i.findtext("z").strip())/lunit
			# get the layer, if not specified, put into default layer
			li=i.findtext("li")
			if li==None:
				li=0
			else:
				li=int(li)
			element=int(i.findtext("el").strip())
			subtype=i.findtext("st")
			if subtype==None:
				subtype=self.PTE[element]
			charge=i.findtext("chr")
			if charge==None:
				charge=0.0
			# now add the atom just read Do not check consistency after each individual atom for performcance reasons
			self.addatom(element, (x,y,z), li, charge, subtype, checkConsistency=False)
		# sanity check now!
		self._consistency_check()




	def readfile(self, filename,typespec=None):
		"""read geometry from specified filename
		known file types are:
			- B{gen} dftb .gen format
			- B{xyz} xmol carthesian format
			- B{fmg} flexible molecular geometry xml format
			- B{cdh} chemical data hierarchy
		@param filename: input file name
		@param typespec: lowercase string, specifying file type, autodetect if omitted
		"""
		typefunctions={
			'gen': self.readgen,
			'xyz': self.readxyz,
			'fmg': self.readfmg,
			'cdh': self.readCDHFile
			}
		if typespec!=None:
			ftype=typespec.strip('.').strip().lower()
		else:
			# strip off compressed format extensions
			if filename.endswith(".gz"):
				tempfilename=filename[:-3]
			elif filename.endswith(".bz2"):
				tempfilename=filename[:-4]
			else:
				tempfilename=filename
			ftype=tempfilename.strip()[-3:].lower()
		if ftype in typefunctions.keys():
			typefunctions[ftype](filename)
		else:
			raise GeometryError("Unknown input file type")



	def read_dftb_charges(self, chargefilename):
		"""reads Mulliken charges from DFTB CHR.DAT file
		@param chargefilename:  input file name
		"""
		if not (os.path.exists(chargefilename) or os.path.exists(chargefilename+".gz") or os.path.exists(chargefilename+".bz2")):
			raise GeometryError("specified charge file does not exist")
		chrfile=utils.compressedopen(chargefilename,'r')
		# skip the first 5 lines
		for i in range(5):
			line=chrfile.readline()  # @UnusedVariable
		atchr=[]
		for i in range(self.Atomcount):
			dummy=chrfile.readline().split()
			#all we want is the atomic charge from the first column
			#we store electron excess as negative charge
			atchr.append(-(float(dummy[1])-self.VALEL[self.AtomTypes[i]]))
		self.AtomCharges=atchr



	def read_noodle_charges(self, detailedoutname):
		"""reads Mulliken charges and l-shell populations from dftb+ detailed.out
		
		I{CAVEAT} Not yet spin-aware!
		@param detailedoutname:  input file name
		"""
		if not (os.path.exists(detailedoutname) or os.path.exists(detailedoutname+".bz2") or os.path.exists(detailedoutname+".gz")):
			raise GeometryError("specified detailed.out file does not exist")
		detfile=utils.compressedopen(detailedoutname,"r")
		detlines=list(detfile)
		detfile.close()
		# find mulliken charges and l-shell population blocks
		chrbase=0
		#lpopbase=0
		line=0
		for i in detlines:
			#atom populations tag
			if i[:18]==" Atom populations ":
				chrbase=line
			#l-shell populations tag
##			elif i[:24]=="  Atom  l         Charge"
##				lpopbase=line
			line+=1
		# incorporate atom charges
		atchr=[]
		# fix for changed output format in DFTB+ 1.1
		# TODO: CHange to reading tagged.out at some future point as detailed.out format is unstable
		if len(detlines[chrbase+1].split())==3:
			# 1.0 format
			lineoffset=1
			populationcolumn=2
		elif len(detlines[chrbase+1].split())==2:
			# 1.1 preliminiary format
			lineoffset=2
			populationcolumn=1
		for i in range(self.Atomcount):
			dummy=detlines[i+chrbase+lineoffset].split()
			#we store electron excess as negative charge
			atchr.append(-(float(dummy[populationcolumn])-self.VALEL[self.AtomTypes[i]]))
		self.AtomCharges=atchr



	def getatomsymlistdict(self):
		"""return a list of atom ordinal numbers present in the geometry and a
		dictionary of ordinal number to atom type number. Both are sorted in 
		order of occurrance"""
		atsymlist=[]
		AtomSymdict=dict()
		for i in self.AtomTypes:
			if i not in AtomSymdict:
				AtomSymdict[i]=len(AtomSymdict)
				atsymlist.append(i)
		return atsymlist,AtomSymdict



	def getFmgString(self):
		"""return a string representing self as fmg geometry ELement"""
		lines=[]
		lines.append('<geometry>')
		lines.append('\t<mode>'+self.Mode+'</mode>')
		if self.Mode=="S":
			lines.append('\t<lattice lunit="ang" orgx="'+str(self.Origin[0]*Angstrom)
				+'" orgy="'+str(self.Origin[1]*Angstrom)
				+'" orgz="'+str(self.Origin[2]*Angstrom)+'">')
			lines.append("\t\t<latvec_a>"+str(self.Lattice[0][0]*Angstrom)+
				" "+str(self.Lattice[0][1]*Angstrom)+" "+str(self.Lattice[0][2]*Angstrom)
				+" </latvec_a>")
			lines.append("\t\t<latvec_b>"+str(self.Lattice[1][0]*Angstrom)+
				" "+str(self.Lattice[1][1]*Angstrom)+" "+str(self.Lattice[1][2]*Angstrom)
				+" </latvec_b>")
			lines.append("\t\t<latvec_c>"+str(self.Lattice[2][0]*Angstrom)+
				" "+str(self.Lattice[2][1]*Angstrom)+" "+str(self.Lattice[2][2]*Angstrom)
				+" </latvec_c>")
			lines.append("\t</lattice>")
		if len(self.LayerDict)>1:
			for i in self.LayerDict.keys():
				lines.append('\t<layer>\n\t\t<li>'+str(i)
					+'</li>\n\t\t<lname>'+self.LayerDict[i].Name
					+'</lname>\n\t</layer>')
		for i in range(self.Atomcount):
			lines.append('\t<atom lunit="ang">\n'+
				"\t\t<x>"+str(self.Geometry[i][0]*Angstrom)+"</x> "
				+"<y>"+str(self.Geometry[i][1]*Angstrom)+"</y> "
				+"<z>"+str(self.Geometry[i][2]*Angstrom)+"</z>\n"
				+"\t\t<el>"+str(self.AtomTypes[i])+"</el>")
			if len(self.LayerDict)>1:
				lines.append("\t\t<li>"+str(self.AtomLayers[i])+"</li>")
			lines.append("\t\t<st>"+self.AtomSubTypes[i]+"</st>")
			lines.append("\t\t<chr>"+str(self.AtomCharges[i])+"</chr>")
			if len(self.LPops[i])>0:
				dummy=[]
				for j in self.Lpops[i]:
					dummy.append("{0: 24.17E}".format(j))
				lines.append("\t\t<lpop>"+"\t".join(dummy)+"</lpop>")
			lines.append("\t</atom>")
		lines.append("</geometry>")
		return "\n".join(lines)

	fmgString=property(getFmgString,doc="String representation of Geometry in .fmg (xml) format")


	def writefmg(self, filename):
		"""write geometry to .fmg file
		@param filename: output file name
		"""
		outstring='<?xml version="1.0" encoding="ISO-8859-1" ?>\n<!DOCTYPE fmg>\n<fmg>\n'
		outstring+=self.fmgString
		outstring+='\n</fmg>\n'
		outfile=open(filename,"w")
		print(outstring,file=outfile)
		outfile.close()



	def writeTurboMole(self,filename="coord"):
		"""write geometry in turbomole coord format
		Does not support definition of internal coordinates!
		@param filename: output filename (default "coord")
		"""
		if self.Mode!="C":
			raise GeometryError("Turbomole file of non-cluster geometry requested")
		outfile=open(filename,"w")
		print(self.tmString(),file=outfile)
		outfile.close()




	def tmString(self):
		"""return a turbomole representation of the geometry
		"""
		#start with coordinates section identifier
		retstr="$coord\n"
		# generate carthesian corrdinales lines using AtomSubType
		for i in range(self.Atomcount):
			# length unit for TM is Bohr
			retstr+="{0: 24.17E} {1: 24.17E} {2: 24.17E} {3:s}\n".format(self.Geometry[i][0]*constants.BOHR,self.Geometry[i][1]*constants.BOHR,self.Geometry[i][2]*constants.BOHR,self.AtomSubTypes[i])
		return retstr



	def getGromacsString(self):
		"""return a string describing self in gromacs file format
		@return: string with geometry representation in gromos87 format"""
		#initialze return string with title and number of atoms
		retstr="comatsci geometry file\n{0:d}\n".format(self.Atomcount,)
		#add atom lines in fixed format as described in gromacs documentation
		for i in range(self.Atomcount):
			# residue name is layer index+1 (to start counting from 1)
			resnum=self.AtomLayers[i]+1
			# take layer name as residue name
			resname=self.LayerDict[self.AtomLayers[i]].Name[0:5].upper()
			# take atom subtype as atom name
			aname=self.AtomSubTypes[i]
			# atom number is index+1
			anum=i+1
			#positions from array, velocities are 0
			retstr+="{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}{7:8.4f}{8:8.4f}{9:8.4f}\n".format(resnum,resname,aname,anum,self.Geometry[i][0]*constants.NANOMETER,self.Geometry[i][1]*constants.NANOMETER,self.Geometry[i][2]*constants.NANOMETER,0.,0.,0.)
		# finally put the cell vectors
		retstr+="{0:12.6f} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f} {5:12.6f} {6:12.6f} {7:12.6f} {8:12.6f}".format(self.Lattice[0][0],self.Lattice[1][1],self.Lattice[2][2],self.Lattice[0][1],self.Lattice[0][2],self.Lattice[1][0],self.Lattice[1][2],self.Lattice[2][0],self.Lattice[2][1])
		#finished, return
		return retstr



	def writegro(self, filename):
		"""write geometry in gromacs/gromos87 .gro file format
		@param filename: output filename
		"""
		outfile=open(filename,"w")
		print(self.getGromacsString(),file=outfile)
		outfile.close()



	def writeAPBS(self, filename):
		"""write geometry data in APBS xml format
		@param filename: name of the output file
		"""
		#construct a root element  "roottag", whoever thought of that element name :-P
		root=ET.Element("roottag")
		# now iterate through atoms and add atom elements...
		for i in range(self.Atomcount):
			currentatom=ET.SubElement(root,"atom")
			# fill current atom with data, lengths are in ANGSTROM
			x=ET.SubElement(currentatom,"x")
			x.text="{0: 24.17E}".format(self.Geometry[i][0]*constants.ANGSTROM,)
			y=ET.SubElement(currentatom,"y")
			y.text="{0: 24.17E}".format(self.Geometry[i][1]*constants.ANGSTROM,)
			z=ET.SubElement(currentatom,"z")
			z.text="{0: 24.17E}".format(self.Geometry[i][2]*constants.ANGSTROM,)
			atomradius=ET.SubElement(currentatom,"atomradius")
			atomradius.text="{0: 24.17E}".format(self.CORAD[self.AtomTypes[i]]*constants.ANGSTROM,)
			charge=ET.SubElement(currentatom,"charge")
			charge.text="{0: 24.17E}".format(self.AtomCharges[i],)
		# now build element tree and write XML to file
		APBSTree=ET.ElementTree(root)
		APBSTree.write(filename)



	def writegen(self, filename, cmode='C'):
		"""Write geometry to dftb .gen file
		possible coordinates modes:
			- B{C} Carthesian coordinates in Angstrom
			- B{F} Fractional coordinates in lattice vector units.
			       Choosing this for a cluster Geometry raises a
			       GeometryError
		@param filename: output file name
		@param cmode: coordinates mode
		"""
		if cmode!='C' and cmode!='F':
			raise GeometryError("unknown coordinate mode")
		if cmode=='F' and self.Mode=='C':
			raise GeometryError("fractional gen for cluster geometry setected")
		outfile=open(filename,'w')
		if cmode=="F":
			writemode="F"
		else:
			writemode=self.Mode
		print(str(self.Atomcount)+"\t"+writemode,file=outfile)
		atlist,AtomSymbols = self.getatomsymlistdict()
		line=""
		for i in atlist:
			line+=self.PTE[i]
			line+=" "
		line+="\n"
		outfile.write(line)
		outfile.flush()
		if self.Mode=="S" and cmode=="F":
			outgeo=self.fractionalGeometry
		else:
			outgeo=self.Geometry*Angstrom
		for i in range(self.Atomcount):
			print("{0:4d} {1:2d}\t{2: 24.17E} {3: 24.17E} {4: 24.17E}".format(i+1, AtomSymbols[ self.AtomTypes[i] ]+1,
				outgeo[i][0], outgeo[i][1], outgeo[i][2]),file=outfile)
		if self.Mode=="S" or cmode=="F":
			print("{0[0]: 24.17E} {0[1]: 24.17E} {0[2]: 24.17E}".format(self.Origin*Angstrom),file=outfile)
			for i in range(3):
				print("{0[0]: 24.17E} {0[1]: 24.17E} {0[2]: 24.17E} ".format(self.Lattice[i]*Angstrom),file=outfile)
		outfile.close()



	def xyzstring(self):
		"""return a string representing self in .xyz format"""
		# to be safe, convert Geometry to an array
		outgeo=numpy.array(self.Geometry,dtype="double")*Angstrom
		# add atom count line
		outstring="{0:d}\n\n".format(self.Atomcount)
		# iterate through atoms and append xyz lines
		for j in range(self.Atomcount):
			outstring+="{0:3}\t{1: 24.17E} {2: 24.17E} {3: 24.17E}\n".format(self.PTE[self.AtomTypes[j]],
			outgeo[j][0],outgeo[j][1],outgeo[j][2])
		return outstring



	def cdfString(self, name, description, atomColumn="e"):
		"""return a string representing self in Balint Aradi's .cdf cluser definitition format
		
		Possible atom columns are:
			- B{e}	element name
			- B{n}	element number
			- B{c}	atomic charge
			- B{s}	subtype (default atom)
		@param name: .cdf name string that uniquely identifies the cluster
		@param description: .cdf description string for the cluster
		@param atomColumn: specifies which information should be listed in the first column of the output file. (default "e")
		@return: string describing self in cluster definition format
		"""
		#sanity check parameters
		if not atomColumn in ("e","n","c","s"):
			raise ValueError("Unknown atomColumn specifier for .cdf format")
		#initialize list of lines to output and insert basic data
		lines=[]
		lines.append("#autogenerated by comatsci")
		lines.append(name)
		lines.append(description)
		lines.append("")
		# output lattice vectors and origin in angstrom
		for i in range(3):
			lines.append(("{24E}"*3).format(tuple(self.Lattice[i]*Angstrom)))
		lines.append("")
		lines.append(("{24E{"*3).format(tuple(self.Origin*Angstrom)))
		lines.append("")
		# output dummy transformation coordinate system
		lines.append("1. 0. 0.")
		lines.append("0. 1. 0.")
		lines.append("0. 0. 1.")
		lines.append("")
		# iterate through coordinates and append atom lines
		for i in range(self.Atomcount):
			# generate first column as specified
			if atomColumn=="e":
				prefix="{0:12s}  ".format(self.PTE[self.AtomTypes[i]])
			elif atomColumn=="n":
				prefix="{0:4d}  ".format(self.AtomTypes[i])
			elif atomColumn=="s":
				prefix="{0:12s}  ".format(self.AtomSubTypes[i])
			elif atomColumn=="c":
				prefix="{0: 24.17E}  ".format(self.AtomCharges[i])
			# generate coordinates in carthesian angstroms
			coordinates=("{24E} "*3).format(tuple(self.Geometry[i]*Angstrom))
			# assemble atom line and append to lines list
			lines.append(prefix+coordinates)
		# finished, convert lines list to output string and return
		return "\n".join(lines)
	
	
	
	def writecdf(self, filename, name, description, atomColumn="e"):
		"""write self to file in Balint Aradi's cluster definition file format
		@see: cdfString for possible atomColumn modes
		@param filename: name of output file to create
		@param name: .cdf name string that uniquely identifies the cluster
		@param description: .cdf description string for the cluster
		@param atomColumn: specifies which information should be listed in the first column of the output file. (default "e")
		"""
		#create output file, build .cdf string of self, write it to output file and close that
		outfile=open(filename,"w")
		print(self.cdfString(name,description,atomColumn),file=outfile)
		outfile.close()




	def writexyz(self, filename, mode='w'):
		"""Writes Geometry to .xyz file. Mode append allows to write animanted .xyz files
		@param filename: 	output file name
		@param mode: file open mode (default w)
		"""
		if mode=='w':
			outfile=open(filename,'w')
		elif mode=='a':
			outfile=utils.compressedopen(filename,'a')
		else:
			raise GeometryError("invalid write mode for .xyz")
		outfile.write(self.xyzstring())
		outfile.close()



	def pointchargesstring(self):
		"""return a string of the current geometry in x y z q format"""
		retstr=""
		outgeo=self.Geometry*Bohr
		for i in range(self.Atomcount):
			retstr+="{0: 24.17E} {1: 24.17E} {2: 24.17E} {3: 24.17E}\n".format(outgeo[i][0],
				outgeo[i][1],outgeo[i][2],self.AtomCharges[i])
		return retstr



	def writexyzq(self, filename):
		"""Write geometry xyzq pointcharges file
		@param filename: output file name
		"""
		outfile=open(filename,"w")
		outfile.write(self.pointchargesstring())
		outfile.close()



	def writepdb(self, filename, occupancy=None, beta=None, writebondlist=0):
		"""Write geometry into pdb file, set ocupancy and beta if specified, otherwise put AtomCharges and zeros respectively. Writes all atoms into HETATM records, writes CONNECT records if writebondlist!=0. <em>CAVEAT: connectivity for periodic structures contains cross-cell bonds!</em>
		@param filename: output file name
		@param occupancy: list of values to put in occupancy field, default: atom charge (default None)
		@param beta: list of values to put in beta files, default: zeroes (default None)
		@param writebondlist: write CONNECT records
		"""
		# open output file
		outfile=open(filename,"w")
		# get PDB string representation of self and print it into the output file
		print(self.getPDBString(occupancy,beta,writebondlist),file=outfile)
		# finished, close output file
		outfile.close()



	def getPDBString(self,occupancy=None, beta=None, writebondlist=0):
		"""return a string containing Geometry in .pdb format. Set ocupancy and beta if specified, otherwise put AtomCharges and zeros respectively. Writes all atoms into HETATM records, writes CONNECT records if writebondlist!=0. <em>CAVEAT: connectivity for periodic structures contains cross-cell bonds!</em>
		@param occupancy: list of values to put in occupancy field, default: atom charge (default None)
		@param beta: list of values to put in beta files, default: zeroes (default None)
		@param writebondlist: write CONNECT records
		"""
		# initialize list of .pdb lines
		pdbLines=[]
		# if no values for occupancy were specified, use atomic occupancies derived from atomic charges
		if occupancy==None:
			occupancy=[]
			for i in range(self.Atomcount):
				occupancy.append(-self.AtomCharges[i]+self.AtomTypes[i])
		elif len(occupancy)!=self.Atomcount:
			raise GeometryError("Occupancy array mismatch")
		# if no beta values were specified, use all zeros
		if beta==None:
			beta=[float(0) for s in range(self.Atomcount)]  # @UnusedVariable
		elif len(beta)!=self.Atomcount:
			raise GeometryError("Beta array mismatch")
		# initialize header
		pdbLines.append("HEADER    comatsci STRUCTURE                        ")
		pdbLines.append("TITLE     comatsci STRUCTURE")
		# check if geometry is supercell and orthorhombic, if so, add CRYST record, else skip supercell!
		cubicmat=numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
		checkmat=self.Lattice*cubicmat
		if self.Mode=="S" and numpy.add.reduce(numpy.add.reduce(checkmat==self.Lattice))==9:
			a=self.Lattice[0][0]*Angstrom
			b=self.Lattice[1][1]*Angstrom
			c=self.Lattice[2][2]*Angstrom
			pdbLines.append("CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}  P 1          1".format(a,b,c,90.,90.,90.))
		# add segment HET records for geometry layers
		for i in self.LayerDict.keys():
			pdbLines.append("HET   {0:3s} {1:4d}".format(self.LayerDict[i].Name[0:3].ljust(3).upper(), i+1))
		# convert coordinates to Angstrom
		outgeo=self.Geometry*Angstrom
		# write HETATM records for each atom
		for i in range(self.Atomcount):
			# ***** begin very very long formatted string *****
			pdbLines.append("HETATM{0:5d} {1:3s}  {2:4s} {3:4d}    {4:8.3f}{5:8.3f}{6:8.3f}{7:6.2f}{8:6.2f}      {9:4s}{10:2s}{11:2s}".format(
				# atom number, Atom Name
				i+1,self.PTE[self.AtomTypes[i]].ljust(3),
				# Layer name (corresponding to HET record created earlier)
				self.LayerDict[self.AtomLayers[i]].Name[0:4].rjust(4).upper(),
				# HET index, X coordinate
				self.AtomLayers[i]+1,outgeo[i][0],
				# Y,Z coordinates
				outgeo[i][1],outgeo[i][2],
				# occupancy and beta as determines above, segment name
				occupancy[i],beta[i],self.LayerDict[self.AtomLayers[i]].Name[0:4].ljust(4).upper(),
				# element name
				self.PTE[self.AtomTypes[i]][0:2].rjust(2).upper(),
				# some empty string
				"  "))
			# *****  end very very long formatted string  *****
		# if bondlist was requested, generate it and append connect records for each bond
		if writebondlist!=0:
			blist=self.bondlist()
			for i in range(self.Atomcount):
				line=copy.deepcopy(blist[i])
				while len(line) > 0:
					token=line[0:4]
					outstring="CONECT{0:5d}".format(i+1)
					for j in token:
						outstring +="{0:5d}".format(j+1)
					pdbLines.append(outstring)
					line=line[4:len(line)]
		# END of geometry entry
		pdbLines.append("END")
		# return lines merged into .pdb format string
		return "\n".join(pdbLines)

	PDBString=property(fget=getPDBString,doc="representation of Geometry in .pdb format")



	def writefdf(self, filename):
		"""Writes Geometry in SIESA .fdf format
		@param filename: output file name
		"""
		outfile=open(filename,'w')
		atlist,AtomSymbols = self.getatomsymlistdict()
		print("NumberOfSpecies {0:d}".format(len(atlist)),file=outfile)
		print("%block ChemicalSpeciesLabel",file=outfile)
		for i in range(len(AtomSymbols)):
			print("{0:5d}\t{1:3d}\t{2:s}".format(i+1,atlist[i],self.PTE[atlist[i]]),file=outfile)
		print("%endblock ChemicalSpeciesLabel\n",file=outfile)
		print("NumberOfAtoms {0:d}\n\nAtomicCoordinatesFormat Ang".format(self.Atomcount),file=outfile)
		print("%block AtomicCoordinatesAndAtomicSpecies",file=outfile)
		outgeo=self.Geometry*Angstrom
		for i in range(self.Atomcount):
			print("{0: 24.17E} {1: 24.17E} {2: 24.17E} {3:d}".format(outgeo[i][0],outgeo[i][1],outgeo[i][2],
					AtomSymbols[self.AtomTypes[i]]+1),file=outfile)
		print("%endblock AtomicCoordinatesAndAtomicSpecies\n",file=outfile)
		if self.Mode=="S":
			print("LatticeConstant  1.0 Ang\n",file=outfile)
			print("%block LatticeVectors",file=outfile)
			outlattice=self.Lattice*Angstrom
			for i in range(3):
				print("{0: 24.17E} {1: 24.17E} {2: 24.17E}".format(outlattice[i][0],outlattice[i][1],
						outlattice[i][2]),file=outfile)
			print("%endblock LatticeVectors",file=outfile)



	def gaussianstring(self):
		"""return geometry specification in carthesian format for Gaussian(03). Does not (yet) write out atom charges."""
		rlines=[]
		tmpgeo=self.Geometry*Angstrom
		for i in range(self.Atomcount):
			rlines.append(" {0:s} {1: 24.17E} {2: 24.17E} {3: 24.17E}".format(self.PTE[self.AtomTypes[i]],
			tmpgeo[i][0],tmpgeo[i][1],tmpgeo[i][2]))
		if self.Mode=="S":
			tmplattice=self.Lattice*Angstrom
			for i in range(3):
				rlines.append(" TV {0: 24.17E} {1: 24.17E} {2: 24.17E}".format(tmplattice[i][0],
					tmplattice[i][1],tmplattice[i][2]))
		return "\n".join(rlines)



	def getHSDChargeConstraintString(self, atomlist, prefactor=1.0):
		"""return a string containing atomic charge constraints in HSD format for DFTB+
		@param atomlist: list of atom indices to provide charge constraints for
		@param prefactor: charge constraint prefactor (default 1.0)
		@return: HSD string specifying charge constraints for the atoms specified in atomlist
		""" 
		# sanity-check atomlist
		for i in atomlist:
			if i >= self.Atomcount:
				raise(ValueError,"atom index out of range")
		# initialize constraints string
		constring=""
		# iterate through atomlist and append the contraint lines
		for i in atomlist:
			constring+="ChargeConstraint = {{ Atoms = {{{0:d}}}\nReferenceCharge = {1: 24.17E}\nPrefactor = {2: 24.17E} }}\n".format(i+1,self.AtomCharges[i],prefactor)
		# finished, return charge constraints string
		return constring
		
		

	def getmasses(self):
		"""return list of atomic masses"""
		m=[]
		for i in range(self.Atomcount):
			m.append(self.AMASS[self.AtomTypes[i]])
		return m



	def compatcheck(self, checkgeo):
		"""Checks, whether given Geometries' atom count, mode, lattice vectors  and type of each atom are compatible
		@param checkgeo: Geometry object to compare
		"""
		if self.Atomcount!=checkgeo.Atomcount:
			raise GeometryError('Atom count mismatch')
		elif self.Mode!=checkgeo.Mode:
			raise GeometryError('Geometry mode mismatch')
		elif not numpy.array_equal(self.Origin, checkgeo.Origin):
			raise GeometryError('Geometry origin mismatch')
		elif not numpy.array_equal(self.Lattice, checkgeo.Lattice):
			raise GeometryError('Geometry lattice mismatch')
		elif self.AtomTypes!=checkgeo.AtomTypes:
			raise GeometryError('Atom type mismatch')



	def setcoordinates(self, newgeometry):
		"""Change coordinates, e.g. for geometry optimization
		@param newgeometry: new coordinates array
		"""
		if numpy.shape(self.Geometry)!=numpy.shape(newgeometry):
			raise GeometryError('Geometry array shape mismatch')
		else:
			self._reset_derived()
			self.Geometry=newgeometry
		self._reset_derived()



	def interpolate(self, other, position):
		"""linerarly interpolate between self and other geomerty at position [position], where position=0 is self and position=1 is other
		@param other: second Geometry for interpolation
		@param position: 0=self, 1=other
		"""
		newgeo=self.__class__(self.Mode, self.Atomcount, self.AtomTypes, self.Origin,
			self.Lattice, self.Geometry, self.AtomLayers, self.LayerDict,
			self.AtomCharges, self.AtomSubTypes)
		if not numpy.array_equal(self.Lattice,other.Lattice):
			newgeo.Lattice=self.Lattice+position*(other.Lattice-self.Lattice)
		if not numpy.array_equal(self.Origin,other.Origin):
			newgeo.Origin=self.Origin+position*(other.Origin-self.Origin)
		newcoords=self.Geometry*(1.0-(float(position)))+other.Geometry*(float(position))
		newgeo.setcoordinates(newcoords)
		return newgeo



	def velcount(self, VALENCE=None):
		"""count valence electrons in geometry
		@param VALENCE: list of valence electron counts, use class defaults if omitted (default None)
		"""
		if VALENCE==None:
			VALENCE=self.VALEL
		elc=0
		for i in self.AtomTypes:
			elc+=self.VALEL[i]
		return elc



	def centerdists(self, center):
		"""return list of the distances of each atom from atom number 'center'
		@param center: origin atom
		"""
		# prepare distances list
		dists=[]
		# sanity checks
		if center < 0 or center > self.Atomcount:
			raise ValueError("central atom index not in atomlist")
		# iterate through atom coordinates
		for i in range(self.Atomcount):
			d=self.Geometry[i]-self.Geometry[center]
			dists.append(numpy.sqrt(numpy.dot(d,d)))
		return dists




	def centerDistSubGeometry(self, center, cutoff):
		"""return a subgeometry object that only contains atoms within cutoff a.u. of center
		@param center: atom index or sequence of three carthesian coordinates marking the center
		@param cutoff: radius of the sphere around center, from which atoms are to be taken
		@return: Geometry object containing only atoms within a sphere of ratius cutoff around center"""
		#first check if an atom index or a sequence of coordinates is supplied as center
		if isinstance(center,(list,tuple,numpy.ndarray)):
			# check coordinates validity
			if len(center)!=3:
				raise ValueError("center coordinates dimmension not equal to 3")
			else:
				centerCoords=numpy.array(center,dtype=float)
		else:
			# check atom index validity
			if not isinstance(center,int):
				raise TypeError("center is not a coordinate sequence or int type")
			if center < 0 or center > self.Atomcount:
				raise ValueError("central atom index not in atomlist")
			else:
				centerCoords=numpy.array(self.Geometry[center],dtype=float)
		# initialize return Geometry object
		returnGeo=self.__class__(self.Mode, iOrigin=self.Origin, iLattice=self.Lattice, iLayerDict=self.LayerDict)
		# calculate distances from center in a new array, then iterate through distance array and append atoms within sphere
		temp=numpy.array(self.Geometry,dtype=float)
		temp-=centerCoords
		distances=numpy.sqrt(numpy.add.reduce(temp*temp,1))
		for i in range(self.Atomcount):
			if distances[i]<=cutoff:
				returnGeo.addatom(self.AtomTypes[i], self.Geometry[i], self.AtomLayers[i], self.AtomCharges[i], self.AtomSubTypes[i], self.LPops[i])
		#finished, return
		return returnGeo




	def bondlist(self, tolerance=1.1):
		"""
		@type tolerance: float
		@param tolerance: tolerance threshold by which to multiply canonical bond lengths when detecting bonds
		@return: list of lists containg bond partners for each atom"""
		if self._blist==None or self._bltolerance!=tolerance:
			self._calcbondlist(tolerance)
			self._bltolerance=tolerance
		return self._blist



	def distancematrix(self):
		"""return the matrix of interatomic distances"""
		if self._dmat==None:
			if self.Mode=="C":
				self._dmat=gx.dmatrix(self.Geometry)
			elif self.Mode=="S":
				self._dmat=gx.sdmatrix(self.Geometry,self.Lattice)
		return self._dmat



	def blmatrix(self):
		"""return the matrix of covalent radius sums"""
		if self._blmat==None:
			#Covalent radii are stored in Angstrom!
			self._blmat=gx.blmatrix(self.AtomTypes,self.CORAD)/Angstrom
		return self._blmat



	def getLinkList(self, othergeo, stepsfunction=None, progressfunction=None):
		"""return a list of links between self and and othergeo
		@param othergeo: Geometry object to find links with
		@param stepsfunction: function to report maximum progress value to, ignore if==None (default None)
		@param progressfunction: function to report progress to, ignore if==None (default None)
		@return: list of tuples of the form (selfindex,otherindex), each representing a bond between self's atom[selfindex] and othergeo's atom[otherindex]
		"""
		#construct a temporary joined geometry
		tempgeo=copy.deepcopy(self)
		tempgeo.appendgeometryflat(othergeo)
		#progress is defined like this here:
		#bondlist calculation=30%, link search=35%, linkatom generation=35%
		if stepsfunction!=None:
			stepsfunction(100)
		#build a list of atom index tuples describing links btw self and othergeo
		bondlist=tempgeo.bondlist()
		if progressfunction!=None:
			progressfunction(30)
		linklist=[]
		for selfatom in range(self.Atomcount):
			for otheratom in bondlist[selfatom]:
				if otheratom >= self.Atomcount:
					linklist.append((selfatom,otheratom-self.Atomcount))
			if progressfunction!=None:
				progressfunction(30+int(35*selfatom/self.Atomcount))
		return linklist







	def _calcbondlist(self,tolerance=1.1):
		"""build list of lists containg bondpartners for each atom,
		calculated from geometry and Covalent radii.
		@type tolerance: float
		@param tolerance: tolearance threshold by which to multiply canonical bond lenths to still detect a bond
		@return: list of atomcount lists, containing bond partner atom indices for each atom"""
		#reimplementation using c extension to calculate all interatomic distances
		#and check for bond on the fly, without having to store and traverse distance
		#matrices. Conserves LOTS of memory for large geometries
		if self.Mode=="C":
			self._blist=gx.blist(numpy.asarray(self.Geometry),self.AtomTypes,self.SBCR, tolerance)
			self._imagecoordlist=None
		elif self.Mode=="S":
			(self._blist,self._imagecoordlist)=gx.sblist(numpy.asarray(self.Geometry), self.Lattice, self.AtomTypes, self.SBCR, tolerance)



	def getSubTypeAtomList(self, subtype):
		"""return a list of atom indices that are of given subtype
		@param subtype: string specifying which atom subtype to extract
		@return: list of atom indices (internal, i.e. counting from zero) that are of _subtype_
		"""
		# initialize list
		stlist=[]
		# iterate through AtomSubTypes list and save indices to return
		for i in range(self.Atomcount):
			if self.AtomSubTypes[i]==subtype:
				stlist.append(i)
		# finished, return
		return stlist



	def getSubTypeSubGeometry(self, subtype, Mode=None):
		"""return a geometry object containing all atoms of self, which are of the specified subtype
		@param subtype: string giving the subtype of atoms to extract
		@param Mode:	mode of the subgeometry object to create. Use
				self.Mode if not specified
				I{Caveat!} Generating a periodic subgeometry of
				a cluster parent will give a useless lattice!
				(default None)
		@return: geometry object with same mode and lattice as self, containing all atoms of _subtype_ in self. Returned geometry onbject will only have the default layer.
		"""
		# if mode is unspecified, use self.Mode
		if Mode==None:
			Mode=self.Mode
		# construct geometry object to return
		stsubgeo=self.__class__(iMode=Mode,iLattice=self.Lattice,iOrigin=self.Origin)
		# get list of subtype atoms and add its atoms to returned subgeometry
		for i in self.getSubTypeAtomList(subtype):
			stsubgeo.addatom(self.AtomTypes[i],self.Geometry[i],None,
				self.AtomCharges[i],self.AtomSubTypes[i])
		# finished, return
		return stsubgeo



	def elemcounts(self):
		"""return a dictionary of element ordinal numbers to numbers of occurrance"""
		symlist,symdict = self.getatomsymlistdict()  # @UnusedVariable
		counts={}
		counts=counts.fromkeys(symlist,0)
		for i in range(self.Atomcount):
			counts[self.AtomTypes[i]]=counts[self.AtomTypes[i]]+1
		return counts



	def distance(self,a,b):
		"""return the distance between atoms a and b
		@param a: first atom index
		@param b: second atom index
		"""
		return self.distancematrix()[a][b]



	def angle(self,a,b,c):
		"""return angle between atoms a-b-c in degrees
		@param a: first atom index
		@param b: apex atom index
		@param c: third atom index
		"""
		# arms of the angle
		arms=numpy.zeros((2,3),dtype=float)
		armlengths=[]
		acidx={0:a,1:c} # helper, indices of a and c to be able to express the periodic expansion as loopoopoops
		for i in (0,1):
				if self.Mode=="S":
						armlengths.append((numpy.sqrt(numpy.dot(self.Lattice[0],self.Lattice[0]))
								+numpy.sqrt(numpy.dot(self.Lattice[1],self.Lattice[1]))
								+numpy.sqrt(numpy.dot(self.Lattice[2],self.Lattice[2]))))
						for u in [-1,0,1]:
								for v in [-1,0,1]:
										for w in [-1,0,1]:
												tempbv=(self.Geometry[acidx[i]]+u*self.Lattice[0]+v*self.Lattice[1]+w*self.Lattice[2])-self.Geometry[b]
												tempbl=numpy.sqrt(numpy.dot(tempbv,tempbv))
												if tempbl<armlengths[i]:
														arms[i]=tempbv
														armlengths[i]=tempbl
		# in a cluster geometry, things are easier...
				else:
						arms[i]=self.Geometry[acidx[i]]-self.Geometry[b]
						armlengths.append(numpy.sqrt(numpy.dot(arms[i],arms[i])))
		# calculate angle from arms and arm lengths
		cosofangle=numpy.dot(arms[0],arms[1])/(armlengths[0]*armlengths[1])
		return math.degrees(math.acos(cosofangle))



	def layerbyname(self,name):
		"""return index of I{first} layer with layer.Name==name, None if
		no match
		@param name: string to match against layer Names
		"""
		for i in self.LayerDict.keys():
			if self.LayerDict[i].Name==name:
				return i
		return None



	def periodicexpand(self,xpns):
		"""periodically expand supercell geometries
		@param xpns: list of 3 integers specifying the new supercell in terms offset original cell vectors"""
		if (self.Mode!="S"):
			raise GeometryError("periodic expansion is only possible for supercells!")
		else:
			# for each lattice vector
			for i in range(3):
				#first allocate an array for the atoms to be appended
				appendlines=numpy.shape(self.Geometry)[0]
				appendgeo=numpy.zeros((appendlines*(xpns[i]-1),3),dtype=float)
				#then calculate the positions of append atoms
				for j in range(1,xpns[i]):
					appendgeo[((j-1)*self.Atomcount):(j*self.Atomcount)]=(
						self.Geometry+(j*self.Lattice[i]))
				#now append the atoms
				for k in range(len(appendgeo)):
					appindex=k%self.Atomcount
					self.addatom(self.AtomTypes[appindex],appendgeo[k],
						self.AtomLayers[appindex],self.AtomCharges[appindex],
						self.AtomSubTypes[appindex])
				#finally update the lattice vector
				self.Lattice[i]=self.Lattice[i]*xpns[i]



	def periodicEmbed(self,shells,surround):
		"""embed self into a shell of periodic geometries
		@param shells: sequence of 3 integers stating how many layers of the surrounding Geometry should be added. <em>This differs from the specification of periodicexpand!</em>
		@param surround: a supercell (mode S) geometry object to apply at the outside of self. If self is mode S, lattices must be identical, otherwise, change self to mode s and apply the lattice of surround
		"""
		#check if surrounding geometry is periodic
		if surround.Mode!="S":
			raise GeometryError("surrounding geometry must be supercell!")
		#if self is not a supercell geometry, make it so, using surrounds lattice
		if (self.Mode!="S"):
			self.Mode="S"
			self.Lattice=numpy.array(surround.Lattice)
			self.Origin=numpy.array(surround.Origin)
		# otherwise, lattices of self and surround must be identical
		elif not numpy.array_equal(self.Lattice,surround.Lattice):
			raise GeometryError("Embedded and surrounding Geometries mut have identical lattices!")
##		# otherwise, self.Lattice vectors must be integer multiples of surround.Lattice
##			else:
##				import numpy.oldnumeric.linear_algebra as LinearAlgebra
##				# calculate the matrix which transforms surround.Lattice into self.Lattice
##				transformMatrix=multiply(self.Lattice,LinearAlgebra.inverse(surround.Lattice))
##				# check if transform matrix is diagonal
##				latticeFactors=diagonal(transformMatrix)
##				if abs(multiply.reduce(latticeFactors)-LinearAlgebra.determinant(transformMatrix))>1E-5:
##					raise ValueError("Incompatible surrounding and core lattices")
##				# check if lattice vector factors are reasonably close to integers
##				for i in range(3):
##					if remainder(latticeFactors[i],1.)>1E-8:
##						raise ValueError("Vector %d of core lattice not integer multiple of surround lattice vector"%i)
##				# if everything is fine, let us 
		# TODO: add check to ensure that embedded and surrounding geometry coordinates do not extend beyond one lattice cell
		# map layers of surround to self: If layer index exists in self, map to that layer, if index does not exist, create new layer:
		for i in surround.LayerDict.keys():
			if not self.LayerDict.has_key(i):
				self.addlayer(surround.LayerDict[i].Name,i)
		# now add the shells of surround
		# loop for each lattice direction a,b,c:
		for a in range(-shells[0],shells[0]+1):
			for b in range (-shells[1],shells[1]+1):
				for c in range (-shells[2],shells[2]+1):
					# self will be in the center, so no need to do anything in that case:
					if a!=0 or b!=0 or c!=0:
						#copy the raw coordinates from surround and shift them
						appendcoordinates=numpy.array(surround.Geometry)
						appendcoordinates+=a*surround.Lattice[0]
						appendcoordinates+=b*surround.Lattice[1]
						appendcoordinates+=c*surround.Lattice[2]
						#now append each atom
						for j in range(surround.Atomcount):
							self.addatom(surround.AtomTypes[j],appendcoordinates[j],surround.AtomLayers[j],surround.AtomCharges[j],surround.AtomSubTypes[j])
		#finally update the lattice vectors
		for i in range(3):
			self.Lattice[i]=self.Lattice[i]*(2.*shells[i]+1.0)



	def totalcharge(self):
		"""return sum of all atomic charges in self"""
		totalcharge=0.0
		for i in range(self.Atomcount):
			totalcharge+=self.AtomCharges[i]
		return totalcharge



	def layerAtomMap(self, layer):
		"""return a dictionary, mapping atom indices inside layer onto global atom indices in self
		@param layer: layer which should be mappend"""
		#the map to return
		lAM={}
		#counter for the intra-layer indices
		intraindex=0
		for i in range(self.Atomcount):
			if self.AtomLayers[i]==layer:
				lAM[i]=intraindex
				intraindex+=1
		return lAM



	def matrixtransform(self, matrix, atomlist=None):
		"""transform positions of the specified atoms by matrix
		@param matrix: the tranformation matrix, no default
		@param atomlist: list of atoms which should be transformed, default: all atoms"""
		if atomlist==None:
			atomlist=range(self.Atomcount)
		for i in atomlist:
			self.Geometry[i]=numpy.dot(matrix,self.Geometry[i])



	def translate(self, tvect, atomlist=None):
		"""translate positions of the specified atoms by tvect
		@param tvect: the vector by which to translate the atom positions, no default
		@param atomlist: list of atoms which should be transformed, default: all atoms"""
		if atomlist==None:
			atomlist=range(self.Atomcount)
		#convert the translation vector to an array
		tarray=numpy.array(tvect,dtype=float)
		for i in atomlist:
			self.Geometry[i]+=tarray




	def foldToCell(self):
		"""Fold coordinates back into periodic unit cell.
		Requires the Geometry instance to have Mode=="S"
		"""
		# check if self is a periodic geometry
		if not self.Mode=="S":
			raise GeometryError("Fold back to unit cell requested on a non-periodic Geometry.")
		# get backfolded coordinates and overwrite current Geometry
		self.setcoordinates(self.getFoldedBackCoordinates())



	def getFoldedBackCoordinates(self):
		"""return array of Geometry coordinates folded back into the unit cell.
		Requires Geometry instance to have mode=="S"
		"""
		# check if self is a periodic geometry
		if not self.Mode=="S":
			raise GeometryError("backfolded coordinates requested on a non-periodic Geometry.")
		# get fractional coordinates
		fractionalCoordinates=self.getFractionalCoordinates()
		# fold back by removing integer components of lattice vectors from fractional coordinates
		# (folds back to fractional coordinates [0...1]. module would fold back to [-1...1])
		foldedFractionalCoordinates=fractionalCoordinates-fractionalCoordinates//1.0
		# project back-folded coordinates back to carthesian space
		foldedCoordinates=numpy.dot(foldedFractionalCoordinates,self.Lattice)
		# finished. return
		return foldedCoordinates



	def getFractionalCoordinates(self):
		"""return array of fractional coordinates of the geometry.
		
		can only be used for Mode=="S" geometries
		"""
		# check if Geometry is in supercell mode
		# GeometryError may be caught and ignored, as unity matrix is used as default Lattice for cluster geometries
		if not self.Mode=="S":
			raise GeometryError("Attempt to calculate fractional coordinates on cluster geometry")
		# calculate inverse of lattice vectors
		inverseLattice=linalg.inv(self.Lattice)
		# project coordintes onto inverse lattice
		fractionalCoordinates=numpy.dot(self.Geometry,inverseLattice)
		# finished, return
		return fractionalCoordinates
	
	fractionalGeometry=property(getFractionalCoordinates,None,None,"Geometry in fractional coordinates. Raises Geometry Error if accessed on cluster Geometries.")
	
	
	
	def centerOnAtom(self, atomIndex):
		""" shift Geometry so that the selected atom lies to the coordinate origin
		@type atomIndex: integer
		@param atomIndex: index of the selected center atom
		"""
		# range check atomIndex
		if atomIndex <0 or atomIndex > self.Atomcount:
			raise ValueError("atom index does not refer to an atom in the geometry")
		else:
			self.translate(self.Geometry[atomIndex]*-1.0, range(self.Atomcount))
		# done
	
	
	
	def rotateAxis(self,axis,angle,atomlist=None,rotatelattice=False):
		""" rotate selected atoms around selected primary axis intersecting coordinate origin
		@type axis: string
		@param axis: single character x,y or z specifying which axis to rotate around
		@type angle: float
		@param angle: rotation angle in raidans
		@type atomlist: sequence of integer
		@param atomlist: indices of atoms to be rotated, if not specified, rotate all atoms
		@type rotatelatiice: boolean
		@param rotatelattlice: if True, also rotate lattice vectors. Useful e.g. to align supercell vectors in x,y,z directions 
		"""
		# check specified axis
		if not axis in ("x","y","z"):
			raise ValueError("invalid rotation axis specified")
		# calculate rotation matrix
		rotmatrix=numpy.zeros((3,3),dtype=float)
		if axis=="x":
			rotmatrix[0][0]=1.0
			rotmatrix[1][1]=math.cos(angle)
			rotmatrix[2][2]=rotmatrix[1][1]
			rotmatrix[2][1]=math.sin(angle)
			rotmatrix[1][2]=-rotmatrix[2][1]
		elif axis=="y":
			rotmatrix[1][1]=1.0
			rotmatrix[0][0]=math.cos(angle)
			rotmatrix[2][2]=rotmatrix[0][0]
			rotmatrix[0][2]=math.sin(angle)
			rotmatrix[2][0]=-rotmatrix[0][2]
		elif axis=="z":
			rotmatrix[2][2]=1.0
			rotmatrix[1][1]=math.cos(angle)
			rotmatrix[0][0]=rotmatrix[1][1]
			rotmatrix[1][0]=math.sin(angle)
			rotmatrix[0][1]=-rotmatrix[1][0]
		# apply rotation matrix to selected atoms
		if atomlist== None:
			self.Geometry=numpy.dot(self.Geometry,rotmatrix)
		else:
			# index-check atomlist
			if (min(atomlist)<0) or (max(atomlist)>self.Atomcount):
				raise ValueError("atom index does not refer to an atom in the geometry")
			for atom in atomlist:
				self.Geometry[atom]=numpy.dot(self.Geometry[atom],rotmatrix)
		if rotatelattice:
			self.Lattice=numpy.dot(self.Lattice,rotmatrix)
		# done
		
	
	def outOfCellPruned(self):
		"""
		return a copy of self with atoms outside the 0th unit cell removed
		@return: Geometry derived from self with any atom outide 0th lattice cell removed
		"""
		pruned=copy.deepcopy(self)
		fractional=pruned.fractionalGeometry
		# construct list of atoms thet need to be removed.
		# list must be in descending order to avoid index shifting problems
		pruneList=[]
		for i in range(pruned.Atomcount-1,-1,-1):
			if (fractional[i][0]>=1.0 or fractional[i][1]>=1.0 or fractional[i][2]>=1.0 or
			 fractional[i][0]<0 or fractional[i][1]<0 or fractional[i][2]<0):
				pruneList.append(i)
		# now remove atoms
		for i in pruneList:
			pruned.delatom(i)
		# finished pruning, return
		return pruned
	
	
	def doubleRemoved(self,doubleCutoff=0.2):
		"""
		return a copy of self with double atoms removed (which may occur at supercell boundaries)
		Only useful for periodic Geometries
		@type doubleCutoff: float
		@param doubleCutoff: if atoms are closer than doubleCutoff * sum of covalent radii, remove one of them  
		@return: Geometry derived from self where atoms falling closer than doubeCutoff to another atom are removed
		"""
		# first get a bond list with short cutoff
		bondList=self.bondlist(doubleCutoff)
		# initialize list of atoms to be removed
		removeList=[]
		# iterate through erroneous short bonds, list one atom to be removed, removed other side of erroneous bond to avoid removing both atoms
		for ii in range(self.Atomcount):
			for jj in bondList[ii]:
				removeList.append(jj)
				bondList[jj].remove(ii)
		# sort list of atoms to be removed
		removeList.sort(reverse=True)
		# make a deep copy of self and remove the offending atoms
		pruned=copy.deepcopy(self)
		for ii in removeList:
			pruned.delatom(ii)
		# finisehd, return
		return pruned
	
		
	
	def layerNames(self):
		"""
		return a list of the layer names in the present geometry
		@return: list of strings
		"""
		return [self.LayerDict[layer].Name for layer in self.LayerDict.keys()]



	@classmethod
	def createMagneli(cls,elementA,elementB,aR,cR,n):
		"""
		create an A_n B_(2n-1) Magneli Phase based on a rutile structure in the axis system defined in
		Le Page and Strobel: J. Sol. Stat. Chem. B{43}, 314 (1982)
		@type elementA: integer
		@param elementA: element number of the cations (metal)
		@type elementB: integer
		@param elementB: element number of the anions (usually oxygen)
		@type aR: float
		@param aR: M{a} lattice constant of the underlying rutile structure
		@type cR: float
		@param cR: M{c} lattice constant of the underlying rutile structure
		@type n: integer
		@param n: M{n} in the AnB2n-1 composition formula. Must be >= 2.
		""" 
		# sanity checks
		if n<2: raise ValueError("n must be >=2")
		# calculate number of B atoms as it will be needed several times
		NB=(2*n)-1
		NBf=float(NB)
		# construct Lattice
		magLattice=numpy.array(
							[[+aR,0.0,-cR],
							 [-aR,-aR,-cR],
							 [0.0,0.0,-((2*n)-1)*cR]])
		# rutile base atoms
		baseElements={'A':elementA,'B':elementB}
		baseCoordinates={'A':numpy.array([[0.,0.,0.],[0.,0.5,0.]]),
		                 'B':numpy.array([[0.0000,0.6952,0.3048],
						                  [0.0000,0.3048,0.6952],
						                  [0.6096,0.8048,0.0856],
						                  [0.3904,0.1952,0.9144]])
						}
		# I-1 symmetry operators
		symmMat=[numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),
				 numpy.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]),
				 numpy.array([[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]]),
				 numpy.array([[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]])
				]
		symmShift=[numpy.array([0.0,0.0,0.0,]),
				   numpy.array([0.5,0.5,0.5,]),
				   numpy.array([0.0,0.0,0.0,]),
				   numpy.array([0.5,0.5,0.5,]),
				]
		# origin shift q (cf. Eq.(2) in Le Page,Strobel)
		if n%2==0:
			q=-0.5
		else:
			q=0.0
		# build list of asymmetric unit cell atoms
		tempCoordinates=[]
		tempTypes=[]
		for k in range(0,(n/2)+1):
			for base in ('A','B'):
				for atomCoord in baseCoordinates[base]:
					Zpk=(atomCoord[2]+(float(k)+q))/NBf
					if Zpk<=0.25 and Zpk>=0.0:
						tempCoordinates.append(numpy.array([atomCoord[0],atomCoord[1],Zpk]))
						tempTypes.append(baseElements[base])
		# apply symmetry operators
		symmCoordinates=[]
		for SO in range(4):
			for coord in tempCoordinates:
				symmCoordinates.append(numpy.dot(coord,symmMat[SO])+symmShift[SO])
		# replicate Types list
		tempTypes=tempTypes*len(symmShift)
		# build geometry instance to return
		returnGeo=cls(iMode="S",iAtomcount=len(tempTypes), iAtomTypes=tempTypes, 
					iLattice=magLattice, iGeometry=numpy.dot(numpy.array(symmCoordinates),magLattice))
		# remove symmetry equivalent atoms before returning
		return returnGeo.doubleRemoved()
		
