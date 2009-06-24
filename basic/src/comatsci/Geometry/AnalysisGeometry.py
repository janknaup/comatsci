## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# AnalysisGeometry.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from Geometry import *
from comatsci import constants,  utils

try:
	from numpy.oldnumeric import linalg
except ImportError:
	from numpy.oldnumeric import linalg

import numpy

import os
import sys
import copy
import math
import bisect


class AnalysisGeometry(Geometry):
	"""Extended Geometry class with geometry analysis functionality"""
	
	# list of features provided by this class
	__geoFeatures__=("analysis", )

	def __init__(self, iMode="C", iAtomcount=0, iAtomTypes=None, iOrigin=None,
				iLattice=None, iGeometry=None):
		"""initialize AnalysisGeometry
		c.f. base class"""
		Geometry.__init__(self, iMode, iAtomcount, iAtomTypes, iOrigin, iLattice, iGeometry)
		self._reset_derived()



	def _reset_derived(self):
		"""reset derived data"""
		Geometry._reset_derived(self)
		self._rbl=None;



	def atomcoordinations(self):
		"""return list of coordination numbers per atom"""
		bl=self.bondlist()
		coordinations=[]
		for i in bl:
			coordinations.append(len(i))
		return coordinations



	def elem_avg_coordinations(self):
		"""return dictionary of average corrdination numbers per element"""
		symlist, symdict = self.getatomsymlistdict()
		elemc={}
		elemc=elemc.fromkeys(symlist,0.0)
		atc=self.atomcoordinations()
		for i in range(self.Atomcount):
			elemc[self.AtomTypes[i]]=elemc[self.AtomTypes[i]]+atc[i]
		counts=self.elemcounts()
		for i in elemc.keys():
			elemc[i]=elemc[i]/counts[i]
		return elemc



	def elem_elem_coordinations(self):
		"""return dictionary of average element-element coordination numbers"""
		elelcoords=self.elem_elem_bondcounts()
		elc=self.elemcounts()
		for i in elc.keys():
			for j in elc.keys():
				elelcoords[i][j]=elelcoords[i][j]/float(elc[i])
		return elelcoords



	def elem_elem_bondcounts(self):
		"""return 2D-dictionary of element-element bond counts"""
		symlist, symdict = self.getatomsymlistdict()
		dummy1={}
		elelc=dummy1.fromkeys(symlist)
		dummy2={}
		for i in elelc.keys():
			elelc[i]=copy.deepcopy(dummy2.fromkeys(symlist,0.0))
		bl=self.bondlist()
		for i in range(self.Atomcount):
			for j in bl[i]:
				elelc[self.AtomTypes[i]][self.AtomTypes[j]]+=1.
		return elelc



	def bondcount(self):
		"""return total number of bonds"""
		bc=0
		bl=self.bondlist()
		for i in bl:
			bc+=len(i)
		return bc/2



	def elem_coordination_counts(self):
		"""return dictionary of dictionaries cotaining counts of coordination numbers by element"""
		symlist, symdict = self.getatomsymlistdict()
		dummy1={}
		el_coord_c=dummy1.fromkeys(symlist)
		for i in el_coord_c.keys():
			el_coord_c[i]={}
		bl=self.bondlist()
		for i in range(self.Atomcount):
			el=self.AtomTypes[i]
			cn=len(bl[i])
			temp=el_coord_c[el].get(cn,0)
			temp += 1
			el_coord_c[el][cn]=temp
		return el_coord_c



	def rt_coordinations(self):
		"""return a string containing coordination statistics in 
		html rich text format"""
		elemavg=self.elem_avg_coordinations()
		elemcounts=self.elemcounts()
		elelco=self.elem_elem_coordinations()
		bc=self.bondcount()
		elcoco=self.elem_coordination_counts()
		lines=[]
		lines.append("<H2>Composition</H2>")
		lines.append("<table border=2 rules=all>\n<tr><th>Element</th><th># of atoms</th><th>%</th></tr>")
		for i in elemcounts.keys():
			workstring="<tr><td>%4s</td><td>%6d</td><td>%5.1f%%</td></tr>" % (self.PTE[i],elemcounts[i],float(elemcounts[i])*100./float(self.Atomcount))
			lines.append(workstring)
		lines.append("<tr><td>%6d</td><td>%5.1f%%</td></tr>\n</table>" % (self.Atomcount,100))
		lines.append("<H2>Coordination</H2>")
		lines.append("<H3>Average coordination by element</H3>")
		lines.append("<table border=2 rules=all>\n<tr><th>Element</th><th>count</th><th>avg. coord.</th></tr>")
		for i in elemavg.keys():
			workstring="<tr><td>%4s</td><td>%6d</td><td>%6.3f</td></tr>" % (self.PTE[i],elemcounts[i],elemavg[i])
			lines.append(workstring)
		lines.append("</table>")
		lines.append("<H3>Element coordination breakdown</H3>")
		lines.append("<table border=2 rules=all>")
		lines.append("<tr><th>Element</th><th>coordination</th><th># of atoms</th><th>% of element</th></tr>")
		for i in elcoco.keys():
			for j in elcoco[i].keys():
				workstring="<tr><td>%2s</td><td>%2d</td><td>%5d</td><td>%5.1f%%</td></tr>"%(self.PTE[i],j,elcoco[i][j],float(elcoco[i][j])/float(elemcounts[i])*100)
				lines.append(workstring)
		lines.append("</table>")
		lines.append("<H3>Average element-element coordinations</H3>")
		lines.append("<table border=2 rules=all>\n<tr><th>Elements</th><th>avg. coord.</th></tr>")
		for i in elelco.keys():
			for j in elelco[i].keys():
				workstring="<tr><td>%2s-%2s</td><td>%6.3f</td></tr>" % (self.PTE[i],self.PTE[j],elelco[i][j])
				lines.append(workstring)
		lines.append("</table>")
		newline="\n"
		return newline.join(lines)



	def reduced_bondlist(self):
		"""return bondlist cleared of double counts"""
		if self._rbl==None:
			self._rbl=copy.deepcopy(self.bondlist())
			for i in range(self.Atomcount):
				for j in (self._rbl[i]):
					self._rbl[j].remove(i)
		return self._rbl



	def rt_bondlist(self):
		"""returns the bondlist as html rich text"""
		bl=self.bondlist()
		bc=self.atom_bondcounts()
		bcols=max(bc)
		lines=['<H2>Bond List</H2>','<table rules=all border=2>']
		lines.append("<tr><th>Atom</th><th colspan=%d> bond partners</th></tr>"%(bcols))
		for i in range(self.Atomcount):
			workstring="<tr><td>%2s%6d</td>" %(self.PTE[self.AtomTypes[i]],i+1)
			for j in range(len(bl[i])):
				at=bl[i][j]
				workstring+="<td>%2s%6d</td>" %(self.PTE[self.AtomTypes[at]],at+1)
			if len(bl[i])<bcols:
				workstring+="<td colspan=%d></td>"%(bcols-len(bl[i]))
			lines.append(workstring+"</tr>")
		lines.append("</table>")
		newline="\n"
		return newline.join(lines)



	def atom_bondcounts(self):
		"""return a list of bond coutns per atom"""
		bl=self.bondlist()
		bc=[]
		for i in bl:
			bc.append(len(i))
		return bc



	def get_atom_coordination_differences(self):
		"""return a list of the differences between atomic coorninations and their standard number of valences (does not account for double bonds)
		@return: per atom array of coordination difference from standard"""
		# first get the atom bond counts and construct an array of standard valence counts
		bondcounts=array(self.atom_bondcounts())
		valences=zeros(self.Atomcount, Float)
		for i in range(len(valences)):
			valences[i]=self.VALENCES[self.AtomTypes[i]]
		# noe return the difference between acutal bond counts and valences
		return bondcounts-valences

	atom_coordination_differences=property(get_atom_coordination_differences, doc="per atom array of coordination difference from standard")



	def histogram(self, data, bins, stepsfunction=None, progressfunction=None):
		"""return a tuple of bins and item counts and outside count, binning data into bins. returnvalue[2] contains the number of items outside the range
		@param data: 1D array or list containing the data to analyze
		@param bins: 1D array or list containing the bins to sort data into. <em>Must be sorted ascending</em>
		@param stepsfunction: callback function to report the total number of progress steps to. For progress display purposes (default None)
		@param progressfunction: callback function to report actual progress to. For progress display purposes. (default None)		
"""
		counts=zeros((len(bins)),Float)
		outside=0
		bincount=len(bins)
		if stepsfunction!=None:
			stepsfunction(len(data))
		progress=0
		for i in data:
			if progressfunction!=None:
				progressfunction(progress)
			#before bisection, check, if datapoint is within binning range
			if i < bins[0]:
				print "low outside"
				outside+=1
			elif i>bins[-1]:
				print "high outside"
				outside+=1
			#find bin using bisection
			else:
				counts[bisect.bisect_left(bins,i)]+=1
			progress+=1
		if progressfunction!=None:
				progressfunction(progress)
		return (bins,counts,outside)



	def elem_avg_charges(self):
		"""return a dictionary of average charges for all elements"""
		symlist, symdict = self.getatomsymlistdict()
		elemchr={}
		elemchr=elemchr.fromkeys(symlist,0.0)
		for i in range(self.Atomcount):
			elemchr[self.AtomTypes[i]]=elemchr[self.AtomTypes[i]]+self.AtomCharges[i]
		counts=self.elemcounts()
		for i in elemchr.keys():
			elemchr[i]=elemchr[i]/counts[i]
		return elemchr



	def elem_charges_hist(self,bincount=10, 
		stepsfunction=None, progressfunction=None, 
		histstepsfunction=None, histprogressfunction=None):
		"""return a dictionary of element charges histograms
		@param bincount: number of bins per histogram
		@param stepsfunction: callback function to report the total number of progress steps to. For progress display purposes (default None)
		@param progressfunction: callback function to report actual progress to. For progress display purposes. (default None)
		@param histstepsfunction: callback function to report the total number of histogram progress steps to. For two-level progress display purposes (default None)
		@param histprogressfunction: callback function to report actual histogram progress to. For two-level progress display purposes. (default None)		
"""
		#first generate an element dictionary of atomic charge lists
		symlist, symdict = self.getatomsymlistdict()
		#this is necessary to get independent empty lists in the dictionary
		#possible python bug?
		elemcl={}.fromkeys(symlist,None)
		for i in elemcl.keys():
			elemcl[i]=[]
		for i in range(self.Atomcount):
			charge=self.AtomCharges[i]
			charr=elemcl[self.AtomTypes[i]]
			charr.append(charge)
			elemcl[self.AtomTypes[i]]=charr
		#give the total number of element histograms as step count
		if stepsfunction!=None:
			stepsfunction(len(elemcl.keys()))
		progress=0
		if progressfunction!=None:
			progressfunction(progress)
		#then bin the atomic charge lists
		for i in elemcl.keys():
			binmin=min(elemcl[i])-1.0e-2
			binmax=max(elemcl[i])+1.0e-2
			binstep=(binmax-binmin)/float(bincount)
			bins=arange(binmin,binmax+binstep,binstep)
			elemcl[i]=self.histogram(elemcl[i],bins,histstepsfunction,histprogressfunction)
			progress+=1
			if progressfunction!=None:
				progressfunction(progress)
		return elemcl



	def rt_charges(self):
		"""return charges statistics on rich text (html) format
		@return: string containing charge statistics html fragment"""
		# Header
		lines=['<H2>Charges</H2>']
		# Total Charges
		lines.append("System total charge: %f e<sup>-</sup>" % (self.totalcharge()))
		# Charges by Layers, if number of layers > 1
		if len(self.LayerDict)>1:
			lines.append("<H3>Charges by Layers</H3>")
			lines.append("<table border=2 rules=all><tr><th>Layer</th><th>charge [e<sup>-</sup>]</th></tr>")
			for i in self.LayerDict.keys():
				tmpgeo=self.layersubgeometry(i)
				lines.append("<tr><td>%s</td><td>%f</td></tr>"%(self.LayerDict[i].Name,tmpgeo.totalcharge()))
			lines.append("</table>")
		# Charges by Elements
		lines.append("<H3>Charge by Elements</H3>")
		lines.append("<table border=2 rules=all><tr><th>Element</th><th>avg. Charge [e<sup>-</sup>]</th></tr>")
		echr=self.elem_avg_charges()
		for i in echr.keys():
			lines.append("<tr><td>%s</td><td>%f</td></tr>"%(self.PTE[i],echr[i]))
		lines.append("</table>")
		# Element-Element charge transfer coefficients
		# get charge transfer coefficients dictionary and create a sorted list of keys
		(dq,S,sigma)=self.getElementElementChargeTransfers()
		dqKeys=dq.keys()
		dqKeys.sort()
		# output as html table
		lines.append("<H3>Charge Transfer Coefficients</H3>")
		lines.append("<table><tr><th>S</th><td>: %12f</td><tr><tr><th>&sigma;<sup>2</sup></th><td>: %12f</td></tr></table>"%(S,sigma))
		lines.append("""<table border="2" rules="all"><tr><th>Elements</th><th>dq<sub>a,b</sub></th></tr>""")
		for i in dqKeys:
			lines.append("<tr><td>%03s-%3s</td><td>%11.8f</td>" % (self.PTE[i[0]],self.PTE[i[1]],dq[i]))
		lines.append("</table>")
		return"\n".join(lines)



	def rdf(self, binwidth=0.2, 
		stepsfunction=None, progressfunction=None, 
		histstepsfunction=None, histprogressfunction=None):
		"""return the radial distribution function over the whole geometry
		@param binwidth: stepwidth for the returned rdf
		@param stepsfunction: callback function to report the total number of progress steps to. For progress display purposes (default None)
		@param progressfunction: callback function to report actual progress to. For progress display purposes. (default None)
		@param histstepsfunction: callback function to report the total number of histogram progress steps to. For two-level progress display purposes (default None)
		@param histprogressfunction: callback function to report actual histogram progress to. For two-level progress display purposes. (default None)
		"""
		#get a 1D array of all bond lengths
		bllist=self.distancematrix().ravel()
		#get the range of the rdf
		rmax=max(bllist)
		#set up the distance bins an get the atom-atom distance histogram
		bins=arange(0, rmax+binwidth,binwidth,Float)
		rdhist=array(self.histogram(bllist,bins,histstepsfunction,histprogressfunction)[1],Float)
		#we will output a 2D array of r,rdf, first put in the distances
		rdf=[bins]
		#now normalize for volume and correct for double counting
##		print rdhist
		#bins are progress, may slow things down a little - style over substance
		if stepsfunction!=None:
			stepsfunction(len(bins-1))
		#for proper rdf normalization calculate the total density
		totaldensity=self.numberDensity
		for i in range(1,len(bins)):
			volume=bins[i]**3-bins[i-1]**3
			volume*=constants.PI*(4.0/3.0)
			#factor 2.0 accounts for double counting
			rdhist[i]/=2.0*volume*totaldensity
			if progressfunction!=None:
				progressfunction(i)
		rdf.append(rdhist)
		#hotfix for inexplicable garbage at r=0: remove r=0 values
		#first remove the bin
		rdf[0]=rdf[0][1:-1]
		#now the values
		rdf[1]=rdf[1][1:-1]
		return rdf



	def getNumberdensity(self):
		"""return the density of the current geometry<br />
		Use supercell volume for periodic structures and and a sphere tightly
		circumcising the molecule otherwise.<em>This breaks for surface models!</em>
		@return:: particle density in the supercell or sphere cuircumscribing the molecule"""
		#refactored volume calculation
		return(float(self.Atomcount)/self.volume)
	
	numberDensity=property(getNumberdensity, doc="particle density in the supercell or sphere cuircumscribing the molecule")



	def getVolume(self):
		"""calculate the volume of the Gemetry, based on supercell for periodic models and circumscribing sphere for clusters
		@return:: volume of the supercell or sphere circumscribing the cluster/molecule"""
		# if self is non-periodic, calculate the volume of the circumscribing sphere
		if self.Mode=="C":
			vol=float(max(self.distancematrix().ravel()))**3*constants.PI*(4.0/3.0)
		else:  #otherwise, calculate the volume of the supercell as (a x b) dot c
			vol=linalg.det(array(self.Lattice))
		# finished, return volume
		return vol

	volume=property(getVolume, doc="volume of the supercell or sphere circumscribing the cluster/molecule")


	def getDipoleMoment(self):
		"""return the total dipole moment of self"""
		# build arrays of coordinates and total charges
		coordinates=array(self.Geometry, Float)
		charges=array(self.AtomCharges, Float)
		# init return vector
		dipole=array((0., 0., 0.,), Float)
		# iterate through arrays (why doesn"t this work on BLAS level?
		for i in range(len(coordinates)):
			dipole+=coordinates[i]*charges[i]
		# finished, return
		return dipole
	
	dipoleMoment=property(getDipoleMoment, doc="Total dipole moment of the geometry's charge distribution")
	
	
	
	def getMassDensity(self):
		"""return the mass density rho of the geometry
		@return:: mass density of self"""
		#calculate total mass
		mass=sum(self.getmasses())
		#return density
		return mass/self.volume
	
	massDensity=property(getMassDensity, doc="mass density of the geometry")
	
	
	
	def getElementElementChargeTransfers(self):
		"""Calculate average charge transfers per bond for each Element-Element combination by least squares method
		@return: (Chargetrans,RMSR) <em>Chargetrans</em>: dictionary of charge transfers, <em>RMSR</em> Root-of-Mean-Squares of Residuals from least squares parametrization.
		"""
		# This method works on the set of equations Q_i=sum_b(dq_ab*c_b), where i runs over all atom indices, a is the i-th atom's element, b runs over all elements present in the geometry, dq_ab is the element-conbination charge transfer factor and c_b is the count of bond partners of type b
		#
		# *** get list of element-element combinations ***
		# first get list of element symbols, sort it and discard unused return values
		symlist,symdict=self.getatomsymlistdict()
		symlist.sort()
		del symdict
		# now construct list of _unique_ combinations and a dictionary mapping element-element pairs to combination list indices
		elementCombinations=[]
		reverseElementCombinations={}
		for i in range(len(symlist)):
			for j in range (i,len(symlist)):
				elementCombinations.append((symlist[i],symlist[j]))
				# put forward and reverse element combinations into reverse dictionary, since both sortings will appear during construction of equations
				reverseElementCombinations[(symlist[i],symlist[j])]=len(elementCombinations)-1
##				reverseElementCombinations[(symlist[j],symlist[i])]=len(elementCombinations)-1
		numElementCombinations=len(elementCombinations)
		# *** construct Q_i equations ***
		# q-factors are sorted like elementCombinations, lines are the equations obtained from each atom
		qMatrix=zeros((self.Atomcount,len(elementCombinations)),Float)
		bondList=self.bondlist()
		for i in range(self.Atomcount):
			iType=self.AtomTypes[i]
			for j in bondList[i]:
				jType=self.AtomTypes[j]
				if iType<jType:
					qMatrix[i][reverseElementCombinations[(iType,jType)]]+=1.
				else:
					qMatrix[i][reverseElementCombinations[(jType,iType)]]-=1.
		# *** get least squares solution of equation system ***
		import numpy.oldnumeric.linear_algebra as LinearAlgebra
		charges=array(self.AtomCharges,Float)
		(q_ab,S,rankA,A)=LinearAlgebra.linear_least_squares(qMatrix,charges)
##		print "q_ab : ",q_ab
##		print "rms  : ",rms
##		print "rankA: ",rankA
##		print "A    : ",A
		# calculate sums of squared residuals by dot product formula
		S=dot(charges,charges)
		S-=dot(dot(qMatrix,transpose(q_ab)),charges)
		# calculate variance (sigma**2)
		variance=S/(float(len(qMatrix))-float(len(q_ab)))
		# *** construct output dictionary of charge transfer coefficients and return ***
		outputQ={}
		q_ab_index=0
		for i in range(len(symlist)):
			for j in range(i,len(symlist)):
				outputQ[(symlist[i],symlist[j])]=q_ab[q_ab_index]
				q_ab_index+=1
		return (outputQ,S,variance)




	def getBondAngles(self):
			"""return a dictionary of all (unique) bond angles in self
			@return: dicitonary with 3-tuple key of atom indices and float values of bond angles in degrees"""
			# initialize bond angles dictionary
			bangdict={}
			# we work on the reduced bond list
			rbl=self.bondlist()
			# iterate through all atoms for a
			for a in range(self.Atomcount):
					# iterate through all bond partners of a as possible apex atom b
					for b in rbl[a]:
							# iterate through all bond partners of apex atom
							for c in rbl[b]:
									if a<c: # bond list assures that a!=b and c!=b, so only a,c must be checked
											bangdict[(a,b,c)]=self.angle(a,b,c)
			# finished, return
			return(bangdict)
	
	
	
	
	def getBondAngleHistogram(self,bins=36,normed=False):
			"""return a histogram of all bond angles in the geometry
			@param bins: number of bins in histogram, default=36 (<=5 deg per bin)
			@param normed: if True, the histogram is normed, so that the area under it is 1. default=False
			@return: array containing bond angles histogram"""
			return numpy.histogram(self.getBondAngles().values(),bins,new=True,normed=normed)




	def getBondAnglesByElementsHistograms(self,bins=36,normed=False):
		"""return a dictionary of histograms of all bond angles in the geometry, grouped by central atom element
		@param bins: number of bins in each histogram, default=36 (<=5 deg per bin)
		@param normed: if True, the histogram is normed, so that the area under it is 1. default=False
		@return: dictionary of arrays containing bond angles histograms, key is elemental Z"""
		# get bond angles dictionary
		bangles=self.getBondAngles()
		# initialize work dictionary
		bangbyelement={}
		# iterate through all bond angles
		for i in bangles.keys():
			cntrat=self.AtomTypes[i[1]]
			# if this angle has new central atom type, add the corresponding list to bangbyelement
			if not cntrat in bangbyelement.keys():
				bangbyelement[cntrat]=[]
			# append bond angle to corresponding
			bangbyelement[cntrat].append(bangles[i])
		# generate histograms for elemental bond angle distributions
		histbyelement={}
		ec=self.elemcounts()
		for i in bangbyelement.keys():
			histbyelement[i]=numpy.histogram(bangbyelement[i],bins,new=True,normed=normed)
			# if histograms are normed, weight each one by the composition fraction of its element, to allow comparbility with the total histogram
			if normed:
				histbyelement[i]=(histbyelement[i][0]/(float(self.Atomcount)/float(ec[i])),histbyelement[i][1])
		# finished, return
		return histbyelement

