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
		bl=self.bondlist(tolerance=1.2)
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
		bl=self.bondlist(tolerance=1.2)
		for i in range(self.Atomcount):
			for j in bl[i]:
				elelc[self.AtomTypes[i]][self.AtomTypes[j]]+=1.
		return elelc



	def bondcount(self):
		"""return total number of bonds"""
		bc=0
		bl=self.bondlist(tolerance=1.2)
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
		bl=self.bondlist(tolerance=1.2)
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
		blstats=self.bondLengthStatistics()
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
		lines.append("<table border=2 rules=all>\n<tr><th>Elements</th><th>avg. coord.</th><th>avg. bond length [A]</th></tr>")
		for i in elelco.keys():
			for j in elelco[i].keys():
				workstring="<tr><td>%2s-%2s</td><td>%6.3f</td>" % (self.PTE[i],self.PTE[j],elelco[i][j])
				if i<=j:
					blkey=(i,j)
				else:
					blkey=(j,i)
				if blstats.has_key(blkey):
					workstring+="<td>%8.4f &plusmn; %8.4f</td>" % (blstats[blkey]["mean"]*constants.ANGSTROM,blstats[blkey]["delta"]*constants.ANGSTROM)
				else:
					workstring+="<td>N/A</td>"
				workstring+="</tr>"
				lines.append(workstring)
		lines.append("</table>")
		newline="\n"
		return newline.join(lines)



	def reduced_bondlist(self,tolerance=1.2):
		"""return bondlist cleared of double counts"""
		if self._rbl==None:
			self._rbl=copy.deepcopy(self.bondlist(tolerance))
			for i in range(self.Atomcount):
				for j in (self._rbl[i]):
					self._rbl[j].remove(i)
		return self._rbl



	def rt_bondlist(self):
		"""returns the bondlist as html rich text"""
		bl=self.bondlist(tolerance=1.2)
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



	def atom_bondcounts(self,tolerance=1.2):
		"""
		@type tolerance: float
		@param toleracne: tolerance threshold by which to multiply canonical bond lengths when detecting bonds
		@return: a list of bond coutns per atom"""
		bl=self.bondlist(tolerance)
		bc=[]
		for i in bl:
			bc.append(len(i))
		return bc



	def get_atom_coordination_differences(self, bondtolerance=1.2):
		"""return a list of the differences between atomic coorninations and their standard number of valences (does not account for double bonds)
		@type bondtolerance: float
		@param bondtolerace: factor applied to canonical bond lengths when counting neighbors as bonded
		@return: per atom array of coordination difference from standard"""
		# first get the atom bond counts and construct an array of standard valence counts
		bondcounts=array(self.atom_bondcounts(bondtolerance))
		valences=zeros(self.Atomcount, Float)
		for i in range(len(valences)):
			valences[i]=self.VALENCES[self.AtomTypes[i]]
		# noe return the difference between acutal bond counts and valences
		return bondcounts-valences

	atom_coordination_differences=property(get_atom_coordination_differences, doc="per atom array of coordination difference from standard")



	def coordinationAnalysis(self, bondtolerance=1.1):
		"""summarize coordination differences by element type
		@type bondtolerance: float
		@param bondtolerance: factor applied to canonical bond lengths when counting neighbors as bonded
		@returntype: dictionary
		@return: dictionary of per-element coordination difference histograms, stored as dictionaries
		
		The format of the returned dictionary is:
		
			- B{stats}
				- B{devMax}: maximum coordination deviation over all elements
				- B{devMin}: minimum (i.e. largest negative) coordination deviation over all elements
			- I{[element]}:
				- I{[deviation]}: number of atoms deviationg by I{deviation} bonds from their canonical valences count
			
		Note that the dictionary is sparse, i.e. no I{deviation} entries for deviation count zero are present.	
		"""
		#get coordination Differences vector
		coordDiffs=self.get_atom_coordination_differences(bondtolerance)
		#initialize return dictionary and statisticsw variables
		analysisDict={}
		for i in self.getatomsymlistdict()[0]:
			analysisDict[i]={}
		deviationMinimum=0
		deviationMaximum=0
		#iterate through coordination difference vector
		for i in range(self.Atomcount):
			if coordDiffs[i]<deviationMinimum: deviationMinimum=coordDiffs[i]
			if coordDiffs[i]>deviationMaximum: deviationMaximum=coordDiffs[i]
			analysisDict[self.AtomTypes[i]][coordDiffs[i]]=analysisDict[self.AtomTypes[i]].get(coordDiffs[i],0)+1
		# append general data to anaysis results
		analysisDict["stats"]={
			"devMax": deviationMaximum,
			"devMin": deviationMinimum
		}
		# finished. return
		return analysisDict




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
		@return: (2,x) array of r,rdf(r)
		"""
		#get a 1D array of all bond lengths
		bllist=self.distancematrix().ravel()
		#****************************************************************
		# UPDATED CODE - USE numpy.histogram instead of homebew function
		# does not provide progress callbacks but is orders of magnitude
		# faster
		#****************************************************************
		# calculate number of bins from bin width and max range
		bincount=int(ceil((max(bllist)+binwidth)/binwidth))
		# pseudo progress function calls before and after histogram binning to
		# stay backwards compatible
		if histstepsfunction!=None:
			histstepsfunction(bincount)
		# generate histogram
		(counts,bins)=numpy.histogram(bllist,bincount,new=True,normed=False)
		if histprogressfunction!=None:
			histprogressfunction(bincount)
#-------------------------------------------------------------------------------
# TODO: rmove old code snippet
#		#get the range of the rdf
#		rmax=max(bllist)
#		#set up the distance bins an get the atom-atom distance histogram
#		bins=arange(0, rmax+binwidth,binwidth,Float)
#		rdhist=array(self.histogram(bllist,bins,histstepsfunction,histprogressfunction)[1],Float)
#		#we will output a 2D array of r,rdf, first put in the distances
#		rdf=[bins]
#		#now normalize for volume and correct for double counting
###		print rdhist
#		#bins are progress, may slow things down a little - style over substance
#-------------------------------------------------------------------------------
		rdf=zeros((2,bincount-1),Float)
		if stepsfunction!=None:
			stepsfunction(bincount-2)
		#for proper rdf normalization calculate the total density
#		totaldensity=self.numberDensity
		#bond length list contains 0 self distances so skip innermost bin
		for i in range(1,len(counts)-1):
			volume=bins[i+1]**3-bins[i]**3
			volume*=constants.PI*(4.0/3.0)
			#factor 2.0 accounts for double counting in distance matrix
#			rdf[1][i]=counts[i]/(2.0*volume*totaldensity)
			rdf[1][i]=counts[i]/(2.0*volume)
			if progressfunction!=None:
				progressfunction(i)
		#bond length list contains 0 self distances so remove r=0 values
		# also, shift r values to fall within center of range bin
		rdf[0]=(bins[1:-1])+binwidth/2.0
		#finished, return
		return rdf



	def elementRDFs(self,binwidth=0.2,
			stepsfunction=None, progressfunction=None, 
			histstepsfunction=None, histprogressfunction=None):
		"""return a dictionary of elemental radial distribution functions
		@param binwidth: stepwidth for the returned rdfs
		@param stepsfunction: callback function to report the total number of progress steps to. For progress display purposes (default None)
		@param progressfunction: callback function to report actual progress to. For progress display purposes. (default None)
		@param histstepsfunction: callback function to report the total number of histogram progress steps to. For two-level progress display purposes (default None)
		@param histprogressfunction: callback function to report actual histogram progress to. For two-level progress display purposes. (default None)		@return: dictionary with keys Z and values elemental rdf (2,x) arrays
		"""
		# get list of elements
		Z,dummy=self.getatomsymlistdict()
		# initialize return dictionary
		elementRDF={}
		for i in Z:
			elementRDF[i]=self.elementsubgeometry(i).rdf(binwidth,
						stepsfunction, progressfunction, 
						histstepsfunction, histprogressfunction)
		# finished, return
		return elementRDF



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




	def getBondAngleStats(self):
		"""compile per central atom element bond angle statistics
		@return dictionary with keys Z and values dictionaries of statistical data
		"""
		# get bond angles
		bondAngles=self.getBondAngles()
		# sort atoms into lists by central element
		anglesByElement={}
		for i in bondAngles.keys():
			centralElement=self.AtomTypes[i[1]]
			if not anglesByElement.has_key(centralElement):
				anglesByElement[centralElement]=[]
			anglesByElement[centralElement].append(bondAngles[i])
		# initialize output dictionary
		outputStats={}
		# calculate statistics for each central element and store in output dictionary
		for i in anglesByElement.keys():
			temparray=numpy.array(anglesByElement[i])
			outputStats[i]={}
			outputStats[i]["mean"]=numpy.mean(temparray)
			sigma=numpy.std(temparray)
			outputStats[i]["sigma"]=sigma
			outputStats[i]["delta"]=sigma/numpy.sqrt(float(len(temparray)))*2.0
		# finished, return.
		return outputStats



	def getElementElementBondlengths(self):
		"""Compile a dictionary of all Element-Element bond lengths
		@return: dictionary with keys (Z,Z) tuples and values lists of bond lengths
		"""
		# get the geometries' reduced bond list (without double counting)
		bondlist=self.reduced_bondlist()
		# imitialize output dictionary
		EEBondLengths={}
		# iterate through all bonds
		for i in range(len(bondlist)):
			for j in bondlist[i]:
				# generate key for this bond (Z1,Z2) with Z1 <= Z2
				if self.AtomTypes[i] <= self.AtomTypes[j]:
					key=(self.AtomTypes[i],self.AtomTypes[j])
				else:
					key=(self.AtomTypes[j],self.AtomTypes[i])
				# initialize list of bond lenghts for (Z1,Z2) if necessary
				if not EEBondLengths.has_key(key):
					EEBondLengths[key]=[]
				# calculate bond length from Geometry (avoid distance matrix)
				# account for bond across periodic boundaries
				if self.Mode=="S":
						length=((sqrt(dot(self.Lattice[0],self.Lattice[0]))
								+sqrt(dot(self.Lattice[1],self.Lattice[1]))
								+sqrt(dot(self.Lattice[2],self.Lattice[2]))))
						for u in [-1,0,1]:
								for v in [-1,0,1]:
										for w in [-1,0,1]:
												tempbv=(self.Geometry[i]+u*self.Lattice[0]+v*self.Lattice[1]+w*self.Lattice[2])-self.Geometry[j]
												tempbl=numpy.sqrt(numpy.dot(tempbv,tempbv))
												if tempbl<length:
														length=tempbl
				# in a cluster geometry, things are easier...
				else:
					d=self.Geometry[i]-self.Geometry[j]
					length=numpy.sqrt(numpy.dot(d,d))
				# append bond length to proper elemental list
				EEBondLengths[key].append(length)
		# finished, return
		return EEBondLengths
	



	def getBondLengthHistograms(self,bins=30):
		"""calculate total and elemental histograms of the bond lengths
		@param bins: Number of bins or array of bin boundaries. Default 30 autoscaled bins
		@return: dictionary of bond length histograms with keys (Z1,Z2) and values arrays of bin centers and bond counts
		special key(-1,-1) contains the total histogram
		"""
		# get the elemental bond length lists and compile a total list of bond lengths
		total=[] 
		elemBondLengths=self.getElementElementBondlengths()
		for i in elemBondLengths.keys():
			total+=elemBondLengths[i]
		# initizalize return dictionary
		bondLengthHistograms={}
		# calculate total bond length histogram
		tempHist=numpy.histogram(total,bins,new=True,normed=False)
		bincenters=tempHist[1][:-1]+((tempHist[1][1]-tempHist[1][0])/2.0)
		bondLengthHistograms[(-1,-1)]=numpy.array([bincenters,tempHist[0]])
		# calculate elemental bond length histograms, using bins from total histogram
		for i in elemBondLengths.keys():
			bondLengthHistograms[i]=numpy.array([bincenters,numpy.histogram(elemBondLengths[i],tempHist[1],new=True,normed=False)[0]])
		# finished, return
		return bondLengthHistograms




	def bondLengthStatistics(self):
		"""calculate statistics on Element-Element bond length distributions
		@return: dictionary with keys (Z1,Z2) element-element combinations and values dictionaries with sting keys and bond length statistical data.
	
		Currently implemented statistics are:
			* B{"mean"} mean bond length
			* B{"sigma"} standard deviation
			* B{"delta"} two standard deviations confidence range
		"""
		# get the elemental bond length lists
		elemBondLengths=self.getElementElementBondlengths()
		# initialize return dictionary
		blstats={}
		# iterate through element-element combinations
		for i in elemBondLengths.keys():
			#initialize statistics dictionary
			blstats[i]={}
			#convert bond length list into array, calculate statistical values and store
			temparray=numpy.array(elemBondLengths[i])
			blstats[i]["mean"]=numpy.mean(temparray)
			sigma=numpy.std(temparray)
			blstats[i]["sigma"]=sigma
			blstats[i]["delta"]=sigma/numpy.sqrt(float(len(temparray)))*2.0
		# finished. return
		return blstats
	



	def RMSD(self,othergeo,**kwargs):
		"""calculate RMS displacement between self and othergeo. If specified, only regard atoms from atomlist
		@return: RMS displacement of specified atoms in atomic units
		@rtype: float
		
		@type othergeo: Geometry
		@param othergeo: Geometry object to calculate RMSD od self against
		@param atomlist: sequence of atom indices (counting from 0) to include in RMSD calculation
		@type atomlist: sequence of integers
		"""
		# check both geometries for compatibility
		self.compatcheck(othergeo)
		# sum square displacements per atom
		alist=kwargs.get("atomlist",range(self.Atomcount)) #: list of atoms to iterate over
		MSD=0.0 #: mean square displacement
		for i in alist:
			dispVect=numpy.array(self.Geometry[i])-numpy.array(othergeo.Geometry[i]) #: displacement vector
			MSD+=numpy.dot(dispVect,dispVect)
		# calculate mean of square displacements and 
		return numpy.sqrt(MSD/float(len(alist)))




	def locateVacancies(self, specvalences=None, tolerance=1.2, ignoreUnknownCanonical=True):
		"""find lattice vacancies by grouping undercoordinated atoms
		
		CAVEAT: CURRENT VERSION IS QUITE DUMB AND IS PRONE TO FAIL IF SEVERAL VACANCIES ARE CLOSE TO EACH OTHER!
		
		TODO: ADD k-means clustering of centroids
		
		@rtype: Geometry subclass of same type as instance that locateVacancies was invoked on
		@return: Geometry object, containing one Atom of Type X for each located vacancy
		
		@type  specvalences: dictionary of integers
		@param specvalences: mapping of element ordinal (Z) to canonical valence count. Atoms with less bond partners will be considered neighbors to a vacancy. Elements, for which no canoncial valence count is given, will use default values, which are probably not useful for transition metals.
		@type  tolerance: float
		@param tolerance: bond length tolerance factor for bond detection. Default 1.2
		@type  ignoreUnknownCanonical: boolean
		@param ignoreUnknownCanonical: if True, ignore atoms that have no known canonical valence count. Default: True
		"""
		# if no special valences are given, construct empty dictionary
		if specvalences==None:
			specvalences={}
		# initialize geometry object to return
		returnGeo=self.__class__(iMode=self.Mode,iLattice=self.Lattice)
		# get bond list
		blist=self.bondlist(tolerance=tolerance)
		#iterate through bond list and add all undercoordinated atoms to undercoordinated atom list
		underList=[]
		for i in range(self.Atomcount):
			# get canonical valence count for current atom, if requested, ignore unknown canonical counts but at least warn.
			if self.AtomTypes[i]<len(self.VALENCES):
				canonicalValence=self.VALENCES[self.AtomTypes[i]]
			else:
				canonicalValence=None
			canonicalValence=specvalences.get(self.AtomTypes[i],canonicalValence)
			if canonicalValence==None:
				if ignoreUnknownCanonical:
					print >> sys.stderr, "WARNING: no canonical valence count for element %s. Ignoring atom Number %d." %(self.PTE[self.AtomTypes[i]],i+1)
					canonicalValences=-1
				else:
					print >> sys.stderr, "ERROR: no canonical valence count for element %s (atom %d). Aborting."%(self.PTE[self.AtomTypes[i]],i+1)
					raise ValueError("no canonical valence found for atom")
			# check if current atom is undercoordinated and if yes, append its index to the list of undercorrdinated atoms
			if len(blist[i])<canonicalValence:
				underList.append(i)
		# if no undercoordinated atoms are found, no further actions are needed. Immediately return empty geometry instance in that case
		if len(underList)==0:
			return returnGeo
		# construct a subgeometry containing only the undercoordinated atoms and calculate their distance matrix
		tempGeo=self.__class__(iMode=self.Mode,iLattice=self.Lattice)
		for atom in underList:
			tempGeo.addatom(self.AtomTypes[atom], self.Geometry[atom], None, self.AtomCharges[atom], self.AtomSubTypes[atom],LPop=None,checkConsistency=True)
		
		print blist
		print underList
		print tempGeo.AtomTypes
		# 
		# Build a list of neighbor atoms belonging to a vacancy. This will work
		# for distant vacancies and can probably serve as a good initial guess
		# for k-means assignment, if that should ever get implemented
		#
		# initialize a boolean mask to mark atoms already assigned to a vacancy
		distMatrix=tempGeo.distancematrix()
		blist=tempGeo.bondlist(tolerance=2.05)
		imagecoordinates=list(tempGeo._imagecoordlist)
		print blist
		print imagecoordinates
		mask=[]
		for i in range(tempGeo.Atomcount): mask.append(True)
		# initialize list of vacancies
		vacancies=[]
		# iterate through upper triangle of distance matrix, assigning atoms to vacancies
		for i in range(tempGeo.Atomcount):     # iterate through all lines of distance matrix
			if mask[i]==True:				   # if atom is already assigned a vacancy, the whole line can be skipped
				mask[i]=False				   # if not, start a new vacancy and mask out atom i
				vacancy=[i]
				for j in range(i+1,tempGeo.Atomcount): # now check all column atoms for proximity (leaving out masked atoms)
					if j in blist[i]:
					##if mask[j] and distMatrix[i,j]< 2.0*( self.SBCR[self.AtomTypes[i]]+self.SBCR[self.AtomTypes[j]])/constants.ANGSTROM:
						mask[j]=False					# add proximate atom and mask it out
						vacancy.append(j)
				vacancies.append(vacancy)			# this vacancy is finished, store it
		print vacancies
		latticeArray=array(tempGeo.Lattice)
		# for each vacancy, calculate the centroid
		for i in range(len(vacancies)):
			# warn about vacancies consisting of exactly two atoms, they may be higher-order bonds
			if len(vacancies[i])==2:
				print >> sys.stderr, "WARNING: Vacancy with two neighbors found, this may actually be a double- or triple bond."
			# skip vacancies consisting of only one atom but warn
			if len(vacancies[i])==1:
				print >> sys.stderr, "WARNING: Isolated undercoordinated atom found. Check Vacancy radius parameter."
			# calculate centroid of vacancy:
			centroid=zeros(3,Float)
			# treat cluster and supercell geometries differently
			if self.Mode=="S":
				# in the supercell, use the vacancy atom's bond partners at their appropriate image coordinates
				coordCount=0
				for atom in vacancies[i]:
					for partner in range(len(blist[atom])):
						coordArray=array(imagecoordinates[atom][partner])
						addLattice=(latticeArray.transpose()*coordArray)
						addvector=add.reduce(addLattice)
						coordCount+=1
						for ii in range(3):
							centroid[ii]+=tempGeo.Geometry[blist[atom][partner]][ii]
							centroid[ii]+=addvector[ii]
				centroid/=coordCount
			else:
				for atom in vacancies[i]:
					for ii in range(3):
						centroid[ii]+=tempGeo.Geometry[atom][ii]
				centroid/=len(vacancies[i])
			returnGeo.addatom(0,centroid,None,0,"X")
		
		returnGeo.foldToCell()
		
		return returnGeo
