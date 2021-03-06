## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# ReactionPath.py
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

import h5py
import numpy as num  # @UnusedImport

from .. import geometry

from ..calculators.calcerror import CalcError  # @UnusedImport
from .. import constants,spline,utils
from .pathiterator import PathIterator, EnergyAccessor, ForcesAccessor


#complicated import statement to make it work with python 2.4 and 2.5
#  see if python 2.5's elementTree implementation is present
try:
	from xml.etree import ElementTree as ET #@UnusedImport
#  otherwise try to import locally installed elementtree (for python 2.4 and below)
except:
	from elementtree import ElementTree as ET #@UnresolvedImport @Reimport


class Reactionpath:
	"""Class for work with reaction paths"""

	def __init__(self,icheckpointdir='checkpoint',ifixedatoms=[],icmode='d',ifmax=1e-4,ifrms=1e-5,imaxit=500,
				charge=0.0,verbosity=None):
		"""initialize Reactionpath
		@param icheckpointdir: directory name for checkpoint storage
		@param ifixedatoms: list of atom indices to keep fixed
		@param icmode: energies and forces calculation mode
			* B{d} external scheduler
		@note: scheduler modes B{s} and B{p} have been removed, paramter icmode is left for backward compatibility
		@param ifmax: max normal force convergence criterion
		@param ifrms: rms normal force convergence criterion
		@param imaxit: maximum number of iterations to perform
		@param charge: system charge in electrons
		@param verbosity: Verbosity level. Choose NORMAL as default. (default None)
		"""
		self.realforces=ForcesAccessor(self)
		self.geos=[]
		self.energies=EnergyAccessor(self)
		self.Atomcount=0
		self.nstep=0
		self.checkpointdir=icheckpointdir
		if verbosity==None:
			self.verbosity=constants.VBL_NORMAL
		else:
			self.verbosity=verbosity
		if ifixedatoms!=None:
			self.fixedatoms=ifixedatoms
		else:
			self.fixedatoms=[]
		if icmode=='d':
			self.realforcesfunc=self.schedcalcenergiesforces
			if self.verbosity >= constants.VBL_NORMAL:
				print("Externally scheduled calculation of energies and forces")
		else:
			raise ValueError("Unknown energies/forces calculation mode")
		self.forcetol=ifmax
		if self.verbosity >= constants.VBL_NORMAL:
			print("ReactionPath: Maximum normal force criterion:      {0:12.6f} a.u.".format(self.forcetol))
		self.rmstol=ifrms
		if self.verbosity >= constants.VBL_NORMAL:
			print("ReactionPath: Maximum RMS normal force criterion:  {0:12.6f} a.u.".format(self.rmstol))
		self.maxit=imaxit
		if self.verbosity >= constants.VBL_NORMAL:
			print("ReactionPath: Maximum number of iterations:        {0:12d}".format(self.maxit))
		self.stopsignal=None
		self.charge=charge
		# spline representation interpolating coordinates, lattice, and charges
		self.splineRep=None
		# Renner subspline representation interpolating only coordinates
		self._rSplineRep=None


	def __iter__(self):
		return PathIterator(self)
	
	
	def __getitem__(self,index):
		if isinstance(index,slice):
			if index.start==None: 
				start=0
			else:
				start=index.start
			if index.step==None:
				step=1
			else:
				step=index.step
			return [ self.__getitem__(ii) for ii in range(start,index.stop,step) ]
		else:
			if index<0: index=index+self.__len__()
			if index < 0 or index >=self.numImages:
				raise IndexError
			else:
				return self.geos[index]

	def __setitem__(self,index,data):
		if isinstance(index,slice):
			if index.start==None: 
				start=0
			else:
				start=index.start
			if index.step==None:
				step=1
			else:
				step=index.step
				dataindex=0
				for ii in range(start,index.stop,step):
					self.geos[ii]=data[dataindex]		
		else:
			if index<0: index=index+self.__len__()
			if index < 0 or index >=self.numImages:
				raise IndexError
			else:
				self.geos[index]=data
				self.splineRep=None


	def getHasSplineRep(self):
		return self.splineRep!=None
	hasSplineRep=property(getHasSplineRep,doc="True if a current spline representation is available")


	def appendgeofile(self, filename,checkCompat=True,TypeSpec=None):
		"""append geometry from input file, try to determine the filetype form filename
		(at this point, only .gen files are supported)
		@param filename: geometry file name
		@param checkCompat: perform compatibility check of all subsequent geometries in the path (default True)
		@param TypeSpec: file type specification string, must be known to Geometry class (default None)
		"""
		temp=geometry.Geometry()
		temp.readfile(filename,TypeSpec)
		self.appendGeoObject(temp,checkCompat)



	def appendGeoObject(self, geo,checkCompat=True):
		"""append geometry object
		@param geo: Geometry Object to Append
		@param checkCompat: perform compatibility check Geometry
		"""
		if self.numimages()!=0 and checkCompat:
			try:
				self.geos[0].compatcheck(geo)
			except geometry.GeometryError,inst:
				if inst.args[0]=='Geometry lattice mismatch':
					print("ReactionPath warning: Geometry lattice mismatch",file=sys.stderr)
				else:
					raise
		else:
			self.Atomcount=geo.Atomcount
			self.amasses=geo.getmasses()
		self.geos.append(geo)
# 		if geo.totalenergy!=None:
# 			if self.has_energies() or self.numimages()==1:
# 				self.energies.append(geo.totalenergy)
		# kill possibly stored spline representation
		self.splineRep=None
		
		
		
	def appendPath(self, appendPath, checkCompat=True):
		"""append a reaction path geometry by geometry
		@type appendPath: Reactionpath instance
		@param appendPath: reactionpath to append to self
		@type checkCompat: boolean
		@param checkCompat: perform compatibility check on appended geometries
		"""
		# iterate through Geometries and append
		for ii in appendPath:
			self.appendGeoObject(ii, checkCompat)
# 		# handle energies and forces
# 		if self.has_energies():
# 			if appendPath.has_energies():
# 				self.energies.extend(appendPath.energies)
# 			else:
# 				print("Reactionpath warning: appending path without energies to path with energies. Dropping existing energies")
# 				self.energies=[]
# 		if self.has_realforces():
# 			if appendPath.has_realsforces():
# 				self.realforces.extend(appendPath.realforces)
# 			else:
# 				print("Reactionpath warning: appending path without forces to path with forces. Dropping existing forces")
# 				self.realforces=[]



	def readXyzPath(self, filename, geoconstructor=geometry.Geometry,progressFunction=None,stepsFunction=None):
		"""read path in xyz format (i.e. A file containing N concatenated xyz geometry strings)
		@param filename: name of the multi-frame xyz file to read
		@param geoconstructor: allows to pass a different geometry object constructor than that of Geometry, e.g. the AnalysisGeometry constructor to analyze a path or trajectory (default Geometry.Geometry)
		"""
		# read the whole xyz file into memory, store as a list and prune trailing empty line
		infile=utils.compressedopen(filename,"r")
		inlist=list(infile)
		infile.close()
		if inlist[-1].strip()=="":
			del inlist[-1]
		# parse atom count from first line and check if total number of lines is consistent with that
		try:
			atomcount=int(inlist[0].strip())
		except:
			print("Atom count '{0:s}' from line 1 of xyz path file '{:s}' could not be parsed. Abort.".format(inlist[1],filename))
			raise
		blocklength=atomcount+2
		if len(inlist)%blocklength!=0:
			raise(ValueError,"Number of lines in input file '{0:s}' does not match atom count. Abort.".format(filename))
		# iterate through xyz blocks, parse and append geometries
		imagecount=(len(inlist))//(atomcount+2)
		if progressFunction != None: progress=progressFunction(total=imagecount)
		for image in range(imagecount):
			temp=geoconstructor()
			try:
				temp.parseXyzString("".join(inlist[image*blocklength:(image+1)*blocklength]))
			except:
				print("Parsing of xyz path image number {0:d} failed. Abort.".format(image+1))
				raise
			try:
				self.appendGeoObject(temp,checkCompat=True)
			except:
				print("Inconsistancy in xyz Path file at image {0:d} detected. Abort.".format(image+1))
				raise
			if stepsFunction != None : stepsFunction(progress,1)
		#finished

#_MANU[
	def readcpmdframes(self, filename,geoconstructor=geometry.Geometry):
		""" this is a copy of readXyzPath combined with the _readforces
		"""
		# read the whole xyz file into memory, store as a list and prune trailing empty line
		infile=utils.compressedopen(filename,"r")
		inlist=list(infile)
		infile.close()
		if inlist[-1].strip()=="":
			del inlist[-1]
		# parse atom count from first line and check if total number of lines is consistent with that
		try:
			atomcount=int(inlist[0].strip())
		except:
			print("Atom count '{0:s}' from line 1 of xyz path file '{0:s}' could not be parsed. Abort.".format(inlist[1],filename))
			raise
		blocklength=atomcount+2
		if len(inlist)%blocklength!=0:
			raise(ValueError,"Number of lines in input file '{0:s}' does not match atom count. Abort.".format(filename))
		# iterate through xyz blocks, parse and append geometries and forces
		imagecount=(len(inlist))//(atomcount+2)
		gradients=[]
		for image in range(imagecount):
			temp=geoconstructor()
			try:
				temp.parseXyzString("".join(inlist[image*blocklength:(image+1)*blocklength]))
			except:
				print("Parsing of xyz path image number {0:d} failed. Abort.".format(image+1))
				raise
			try:
				self.appendGeoObject(temp,checkCompat=True)
			except:
				print("Inconsistancy in xyz Path file at image {0:d} detected. Abort.".format(image+1))
				raise
			gradsbuf=[]
			for line in inlist[(image*blocklength+2):(image+1)*blocklength]:
				buf=line.split()
				gradsbuf.append([ float(s)/0.529177 for s in buf[4:7] ])
			gradients.append(numpy.array(gradsbuf))		
			self.realforces=gradients  #FIXME: needs to ahdere to new logic of storing forces in Geometry
								
		#finished
#_MANU]


	def pathinterpolate(self, nimages):
		"""linearly interpolate nimages intermediate images 
		between every two existing path images
		@param nimages: number images to interpolate between every two existing images
		"""
		#newimagecount=len(self.geos)+(len(self.geos)-1)*nimages
		newgeoarray=[]
		for i in range(len(self.geos)-1):
			newgeoarray.append(self.geos[i])
			for j in range(nimages):
				newgeoarray.append(self.geos[i].interpolate(self.geos[i+1],float(j+1)/float(nimages+1)))
		newgeoarray.append(self.geos[len(self.geos)-1])
		self.geos=newgeoarray



	def numimages(self):
		"""return the number of images in the path"""
		return len(self.geos)
	numImages=property(numimages,doc="""number of images in path""")
	def __len__(self):
		return len(self.geos)


	def writegenpath(self,nameprefix="path",minimaloutput=False):
		"""write the path geometries as gen files
		@param nameprefix: string to prepend to each filename
		"""
		for i in range(self.numimages()):
			filename="{0:s}.{1:04.2f}.gen".format(nameprefix, (float(i)/(self.numimages()-1)))
			self.geos[i].writegen(filename)



	def writexyzpath(self,name="path.xyz",minimaloutput=False):
		"""write the path geometries as a single multiframe .xyz file
		@param name: output filename
		"""
		self.geos[0].writexyz(name)
		for i in range(1,self.numimages()):
			self.geos[i].writexyz(name, mode="a")


	
	def fmgFEstring(self):
		"""return a string containing .fmg <trjstep> elements describing forces and energies along the path"""
		FElines=[]
		for i in range(self.numimages()):
			FElines.append("<trjstep>\n\t<nrg>"+str(self.energies[i])+"</nrg>")
			if len(self.realforces)>0:
				FElines.append("\t<forces>")
				for j in range(self.geos[0].Atomcount):
					FElines.append("\t\t{0: 24.17E}  {1: 24.17E}  {2: 24.17E}".format(self.realforces[i][j][0],
						self.realforces[i][j][1],self.realforces[i][j][2]))
				FElines.append("\t</forces>")
			FElines.append("</trjstep>")
		return "\n".join(FElines)



	def fmgTrjinfoString(self):
		"""return a string containing .frm <trjinfo> element 
		giving additional trajectory info, as long as there is any info to return"""
		infolines=[]
		#only append stepcount if >0
		if self.nstep>0:
			infolines.append("\t<stepcount>{0:d}</stepcount>".format(self.nstep))
		#only return <trjinfo> element, if not empty
		if len(infolines)>0:
			return "<trjinfo>\n"+"\n".join(infolines)+"\n</trjinfo>"
		else:
			return ""



	def writefmgpath(self,name='path.fmg',minimaloutput=False):
		"""write the path as a single, multi-geometry .fmg file
		@param name: output file name
		"""
		outstring='<?xml version="1.0" encoding="ISO-8859-1" ?>\n<!DOCTYPE fmg>\n<fmg>\n'
		for i in range(self.numimages()):
			outstring+=(self.geos[i].fmgString+"\n")
		if self.has_energies():
			outstring+=(self.fmgFEstring()+"\n")
		outstring+=self.fmgTrjinfoString()+"\n"
		outstring+='</fmg>\n'
		outfile=open(name,"w")
		print(outstring,file=outfile)
		outfile.close()

	
	

	def etReadFmgPath(self, filename,checkCompat=True,GeoConstructor=geometry.Geometry):
		"""read the path from a single, multi-geometry .fmg file using ElementTree
		@param filename: 	input file name
		@param checkCompat: perform compatibility check of all subsequent geometries in the path (default True)
		@param GeoConstructor: allows to pass a different geometry object constructor than that of Geometry, e.g. the AnalysisGeometry constructor to analyze a path or trajectory (default Geometry.Geometry)
		"""
		# initialize empty geometries list
		self.geos=[]
		# parse xml file and construct overall ElementTree
		if self.verbosity >= constants.VBL_DEBUG2:
			print('Parse xml file "{0:s}" and build elementTree... '.format(filename),end="")
			sys.stdout.flush()
		# handle the file outside of ET to allow use of transparent decompression
		fmgfile=utils.compressedopen(filename)
		tree = ET.parse(fmgfile)
		fmgfile.close()
		if self.verbosity >= constants.VBL_DEBUG2:
			print('done.')
		# iterate through Geometry Elements
		if self.verbosity >= constants.VBL_DEBUG2:
			print('iterate through geometries ',end="")
		geometryElements=tree.findall("geometry")
		for currentGeometryElement in geometryElements:
			tempGeometry=GeoConstructor()
			tempGeometry.handleGeoElementTree(currentGeometryElement)
			if self.verbosity >= constants.VBL_DEBUG2:
				print('.',end="")
				sys.stdout.flush()
			if self.numimages()>1 and checkCompat:
				try:
					self.geos[0].compatcheck(tempGeometry)
				except geometry.GeometryError,inst:
					if inst.args[0]=='Geometry lattice mismatch' and self.verbosity>=constants.VBL_SILENCE:
						print("ReactionPath warning: Geometry lattice mismatch")
					else:
						raise
			self.geos.append(tempGeometry)
		# geometries list build, now get some path global Data
		self.amasses=self.geos[0].getmasses()	
		self.Atomcount=self.geos[0].Atomcount
		# check if trajectory steps data is present and store if true
		trjsteps=tree.findall("trjstep")
		if len(trjsteps)!=self.numimages() and len(trjsteps)!=0 and self.verbosity>=constants.VBL_SILENCE:
			print("ReactionPath warning: Pathfile contains trajectory steps but number does not match image count.\n ignoring forces and energies from file.")
			print(len(trjsteps))
			print(self.numimages())
		elif len(trjsteps)==self.numimages():
			if self.verbosity >= constants.VBL_DEBUG2:
				print('Interate through trajectory steps ',end="")
			#self.energies=[]
			#self.realforces=[]
			self._stepIndex=0
			for i in trjsteps:
				self.handletrjstep_etree(i)
				if self.verbosity >= constants.VBL_DEBUG2:
					print('.',end="")
				self._stepIndex=self._stepIndex+1
		if self.verbosity >= constants.VBL_DEBUG2:
			print('done.')
		# check if global trahectory info is present
		trjinfo=tree.findall("trjinfo")
		if len(trjinfo)!=0:
			if self.verbosity >= constants.VBL_DEBUG2:
				print('Handle trajectory info.')
			self.handletrjinfo_etree(trjinfo[0])
		# kill possibly stored spline representation
		self.splineRep=None
		if self.verbosity >= constants.VBL_DEBUG2:
			print('done reading .fmg path.')


	readfmgpath=etReadFmgPath


	def handletrjinfo_dom(self, info):
		"""handle trajectory info from .fmg file:
		If spresent, store initial stepcount and increase iteration maximum
		@param info: dom object containing the trjinfo subtree
		"""
		#handle stepcount, ignore if not present
		stepcount=info.getElementsByTagName("stepcount")
		if len(stepcount)==1:
			self.nstep=int(stepcount[0].childNodes[0].data)
			#we must, of course, increase the maximum iteration count by 
			#the initial offset applied to the stepcount
			self.maxit+=self.nstep
		elif len(stepcount)>1:
			raise "FMG ERROR: multiple trajectory stepcounts encountered"




	def handletrjinfo_etree(self, info):
		"""handle trajectory info from .fmg file using ElementTree:
		If spresent, store initial stepcount and increase iteration maximum
		@param info: dom object containing the trjinfo subtree
		"""
		#handle stepcount, ignore if not present
		stepcount=info.find("stepcount")
		if stepcount!=None:
			self.nstep=int(stepcount.text)
			#we must, of course, increase the maximum iteration count by 
			#the initial offset applied to the stepcount
			self.maxit+=self.nstep




	def handletrjstep_dom(self, step):
		"""handle trajectory steps from .fmg file and store energies and forces"""
		#first get the total energy
		energies=step.getElementsByTagName("nrg")
		#no default value here:
		if len(energies)!=1:
			raise "FMG ERROR: trajectory step for path must contain exactly one energy!"
		else:
			#self.energies.append(float(energies[0].childNodes[0].data))
			self.energies[self._stepIndex]=float(energies[0].childNodes[0].data)
		# now get the forces, no defaults here either
		forces=step.getElementsByTagName("forces")
		if len(forces)>1:
			raise "FMG ERROR: trajectory step for path must contain exactly one forces block!"
		elif len(forces)<1:
			print("FMG WARNING: no forces in trajectory step")
		else:
			dummy=forces[0].childNodes[0].data.strip().split()
			gradsbuf=[ float(s) for s in dummy ]
			#self.realforces.append(numpy.reshape(numpy.array(gradsbuf),(-1,3)))
			self.realforces[self._stepIndex]=numpy.reshape(numpy.array(gradsbuf),(-1,3))



	def handletrjstep_etree(self, step):
		"""handle trajectory steps from .fmg file and store energies and forces using ElementTree"""
		#first get the total energy
		energies=step.findall("nrg")
		#no default value here:
		if len(energies)!=1:
			raise "FMG ERROR: trajectory step for path must contain exactly one energy!"
		else:
			#self.energies.append(float(energies[0].text))
			self.energies[self._stepIndex]=float(energies[0].text)
		# now get the forces, no defaults here either
		forces=step.findall("forces")
		if len(forces)>1:
			raise "FMG ERROR: trajectory step for path must contain exactly one forces block!"
		elif len(forces)<1:
			print("FMG WARNING: no forces in trajectory step")
		else:
			dummy=forces[0].text.strip().split()
			gradsbuf=[ float(s) for s in dummy ]
			#self.realforces.append(numpy.reshape(numpy.array(gradsbuf),(-1,3)))
			self.realforces[self._stepIndex]=numpy.reshape(numpy.array(gradsbuf),(-1,3))



	def writeCDHPath(self,filename="path.cdh",savespace=True,progressFunction=None,stepsFunction=None,minimaloutput=False):
		""" 
		Write the current path in HDF5 format according to CDH specification
		@type filename: string
		@param filename: name of the HDF5 file to create
		@type savespace: boolean
		@param savespace: try to save space in .cdh file by moving constant data to globals group
		@type progressFunction: callable
		@param progressFunction: function to call after every write step to indicate progress to tge user
		@type stepsFunction: callable taking integer argument
		@param stepsFunction: function to call before writing, to initialize progress indicator.
			Must accept integer argument of total step count.
		"""
		# open HDF5 file for overwriting
		pathfile=h5py.File(filename,"w")
		if stepsFunction!=None: stepsFunction(self.numimages())
		# iterate through path images
		for image in range(self.numimages()):
			imagelabel="frame{0:010d}".format(image)
			if savespace:
				try:
					globalSets=["elements","types","lattice","residues","method"]
					self.geos[0].compatcheck(self.geos[image])
				except geometry.GeometryError,inst:
					if inst.args[0]=='Geometry lattice mismatch':
						print("ReactionPath warning: Geometry lattice mismatch")
						globalSets.remove("lattice")
					else:
						savespace=False
						globalSets=[]
			# first write the image geometry
			if image==0 and savespace:
				globalExclude=set(self.geos[0].knownCDHFields).difference(set(globalSets))
				globalsGroup=self.geos[image].writeCDHFrameGroup(h5file=pathfile,groupname="globals",exclude=globalExclude)[1]  # @UnusedVariable
			if minimaloutput:
				tempMinimalOutput = ["coordinates","elements","types","lattice"]
				excludeMinimal = list(set(list(self[0].knownCDHFields)).difference(set(tempMinimalOutput)))
			else:
				excludeMinimal = []
			imagegroup=self.geos[image].writeCDHFrameGroup(h5file=pathfile,groupname=imagelabel,exclude=globalSets+excludeMinimal)[1] #@UndefinedVariable @UnusedVariable
			if progressFunction: progressFunction()
		pathfile.close()


	def readCDHPath(self,filename,checkCompat=True,geoconstructor=geometry.Geometry,progressFunction=None,stepsFunction=None):
		"""
		Read path from HDF5 file according to cdh specification
		@type filename: string
		@param filename: name of the HDF5 file to read
		@return: reference to self
		"""
		#open HDF5 file for reading
		pathfile=h5py.File(filename,"r")
		# get list of CDH frames, filter out non-frame data sets and sort by index number in name
		sets=pathfile.keys()
		frames=[]
		for ii in sets:
			if ii[0:5]=="frame":
				frames.append(ii)
		frames.sort()
		# iterate through frames
		if stepsFunction!=None: stepsFunction(len(frames))
		globalsGroup=pathfile.get("globals", None)
		for frame in frames:
			tg=geoconstructor()
			framegroup=pathfile[frame]
			tg.parseH5Framegroup(framegroup,globalsGroup)
			self.appendGeoObject(tg, checkCompat)
			if progressFunction!=None: progressFunction()
		
		pathfile.close()
		return(self)


	
	def has_energies(self):
		"""return true if calculated energies are available, else false"""
		#return len(self.energies)==self.numimages()
		for ii in range(self.numImages):
			if self.energies[ii]==None:
				return False
		return True
	
	
	
	def has_realforces(self):
		"""return true if calculated forces are available, else false"""
# 		if self.realforces!=None:
# 			return len(self.realforces)==self.numimages()
# 		else:
# 			return False
		for ii in range(self.numImages):
			if self.realforces[ii]==None:
				return False
		return True



	def _writeforces(self,name="neb.frc"):
		"""write the path real forces to a file
		@param name: output filename
		"""
		outFile=open(name,"w")
		print("{0:6d} {1:6d}".format(self.Atomcount, self.numimages()),file=outFile)
		for i in range(self.numimages()):
			print("{0:d}".format(i),file=outFile)
			tempforce=self.realforces[i].ravel()
			for j in range(self.Atomcount):
				print("{0: 24.17E} {1: 24.17E} {2: 24.17E} ".format(tempforce[3*j],tempforce[(3*j)+1],tempforce[(3*j)+2]),file=outFile)
		outFile.close()



	def _writeenergies(self,name="neb.nrg"):
		"""write the path energies to a file
		@param name: output filename
		"""
		outFile=open(name,"w")
		for i in range(self.numimages()):
			print("{0:5d}  {1: 24.17E}".format(i,self.energies[i]),file=outFile)
		outFile.close()



	def _readforces(self,name="neb.frc"):
		"""read the forces from .frc dump file into path realforces
		@param name: input filename
		"""
		inFile = utils.compressedopen(name,"r")
		buf=inFile.readline()
		dummy=buf.split()
		if self.Atomcount!=int(dummy[0]):
			raise "Atom count mismatch in forces file"
		if self.numimages()!=int(dummy[1]):
			raise "Image count mismatch in forces file"
		gradients=[]
		for i in range(self.numimages()):
			buf=file.readline()
			if int(buf)!=i:
				raise("Error reading force file")
			gradsbuf=[]
			for j in range(self.Atomcount): #@UnusedVariable
				buf=file.readline()
				dummy=buf.split()
				gradsbuf.append([ float(s) for s in dummy[0:3] ])
			gradients.append(numpy.array(gradsbuf))
		for ii in range(self.numimages()): self.realforces[ii]=gradients[ii]
		file.close()



	def _readenergies(self,name="neb.nrg"):
		"""read energies from .nrg dump file
		@param name: input filename
		"""
		inFile=utils.compressedopen(name,"r")
		enbuf=[]
		for line in inFile:
			dummy=line.split()
			enbuf.append(float(dummy[1]))
		if len(enbuf)!=self.numimages():
			raise "Image count mismatch in energy file"
		else:
			inFile.close()
			for ii in range(len(enbuf)): self.energies[ii]=enbuf[ii]



	def readenergiesforces(self,efile,ffile):
		"""Read Energies and forces from dump files
		@param efile: energies input filename
		@param ffile: forces input filename
		"""
		self._readenergies(efile)
		self._readforces(ffile)



	def writecheckpoint(self, directory):
		"""write a path checkpoint into directory
		@param directory: checkpoint directory name
		"""
		curdir=os.path.realpath(".")
		if not os.path.exists(directory):
			os.mkdir(directory)
		os.chdir(directory)
		self.writegenpath()
		if self.has_energies():
			self._writeenergies()
		if self.has_realforces():
			self._writeforces()
		os.chdir(curdir)
		self.writefmgpath(directory+".fmg")


	
	def readcheckpoint(self, checkpointdir,checkCompat=True):
		"""read geometries, energies and forces from checkpoint directory
		@param checkpointdir: checkpoint directory name
		@param checkCompat: perform compatibility check of all subsequent geometries in the path (default True)
		"""
		if os.path.isdir(checkpointdir):
			if self.verbosity >= constants.VBL_NORMAL:
				print("Reading path images, energies and forces from checkpoint {0:s}".format(checkpointdir))
			imglist=[]
			dirlist=os.listdir(checkpointdir)
			for i in dirlist:
				if (i.endswith('.gen') or i.endswith('.gen.gz') or i.endswith('.gen.bz2'))and i.startswith('path'):
					imglist.append(i)
			imglist.sort()
		else:
			errstr= "Checkpoint directory {0:s} does not exist".format(checkpointdir)
			raise ValueError(errstr)
		for i in imglist:
			self.appendgeofile(checkpointdir+'/'+i,checkCompat=checkCompat)
		if self.verbosity >= constants.VBL_NORMAL:
			print("Sucessfully read {0:d} images into path".format(self.numimages()))
		if checkCompat:
			self.readenergiesforces(checkpointdir+'/'+'neb.nrg',checkpointdir+'/'+'neb.frc')
		else:
			self._readenergies(checkpointdir+'/'+'neb.nrg')
		if self.verbosity >= constants.VBL_NORMAL:
			print("Forces and Energies read successfully")
		# kill possibly stored spline representation
		self.splineRep=None




	def _fixatoms(self,force):
		"""return force with zero forces for atoms in self.fixedatoms
		@param force: input forces array
		"""
		rforce=[]
		fshape=numpy.shape(force[0])
		for i in range(self.numimages()):
			tempforce=numpy.reshape(force[i],(self.Atomcount,3))
			for j in self.fixedatoms:
				tempforce[j-1]=numpy.zeros((3),dtype=float)
			rforce.append(numpy.reshape(tempforce,fshape))
		return rforce



	def importEnergiesForces(self, calculator, directories):
		"""Import results from external calculation runs
		@param calculator: calculator object to use for reading results
		@param directories: sequence strings storing directory paths
		"""
		if not len(directories)==self.numimages():
			raise ValueError,"Number of images does not match number of results paths"
		for i in range(self.numimages()):
			calculator.parseForeignOutput(self.Atomcount,directories[i])
			self.energies[i]=calculator.getenergy()
			self.realforces[i]=calculator.getforces()
			calculator.finreready()
			if self.verbosity>=constants.VBL_NORMAL:
				if i%10==0:
					sys.stdout.write(':')
				else:
					sys.stdout.write(".")
				sys.stdout.flush()



	def calcenergiesforces(self, calculator, charge=0, steplabelprefix="step"):
		"""Use calculator to calculate energies and forces
		@param calculator: calculator object
		@param charge: charge for calculation (default 0)
		@param steplabelprefix: prefix prepend to step label
		"""	
		raise NotImplementedError


	def pcalcenergiesforces(self, scheduler, charge=0, steplabelprefix="step"):
		"""Calculate Energies and forces in parallel mode
		calculators have already been initialized on the slaves
		@param scheduler: parallel calculation scheduler function
		@param charge: charge to pass to calculators
		@param steplabelprefix: string to prepend to step labels
		"""
		raise NotImplementedError



	def schedcalcenergiesforces(self, scheduler, charge=0, steplabelprefix="step"):
		"""Calculate Energies and forces, using comatsci generic scheduler class
		@param scheduler: scheduler instance
		@param charge: charge to pass to calculators
		@param steplabelprefix: string to prepend to step labels
		"""
		if not (self.has_energies() and self.has_realforces()):
			start=0
			end=0
			#First energies & forces call, so we must do all images
			#self.energies, self.realforces=self.__schedperform(charge, steplabelprefix, scheduler,0,0)
		else:
			#Subsequent call, we don't have to calculate for the first and last image again, they don't change
			start=1
			end=1
		tempEnergy, tempForce=self.__schedperform(charge, steplabelprefix, scheduler,start,end)
		for ii in range(start,self.numImages-end):
			self.energies[ii]=tempEnergy[ii-start]
			self.realforces[ii]=tempForce[ii-start]



	def __schedperform(self, charge, steplabelprefix, scheduler,startoff,endoff):
		"""helper function to let the sheduler perform energies and forces calculations
		@param charge: Charge to pass to the calculator
		@param steplabelprefix: prefix for the autogenerated steplabels
		@param scheduler: scheduler instance to use
		@param startoff: offset at start of images array, i.e. the first image to calculate for. Must be 0 or positive.
		@param endoff: offset at end of images array, i.e. number of images to ignore at and of images array, Must be 0 or positive.
		"""
		# assemble options dictionaries
		joblist=[]
		for i in range(startoff,self.numimages()-endoff):
			label="{0:s}-{1:6.4f}".format(steplabelprefix,(float(i)/(float(self.numimages())-1)))
			joblist.append({
				"Geometry":copy.deepcopy(self.geos[i]),
				"Charge":charge,
				"steplabel":label,
		})
		# now let the scheduler perform the jobs on the list
		results=scheduler.perform(joblist)
		# now unpack the results into separate arrays
		energies=[]
		forces=[]
		for i in results:
			energies.append(i["Energy"])
			forces.append(i["Forces"])
		# return the results
		return energies, forces



	def rmsforce(self,images=None,force=None):
		"""Calculate the rms force for images, use whole path if omitted
		@param images: if specified, only regard images in this list, otherwise use all images
		@param force: if specified, calculate rms of this forces array, otherwise use NEB forces
		"""
		if force==None:
	#TODO rmsforce -> realforces
			force=self.realforces
		ms=0.0
		if images==None:
			rng=range(self.numimages())
		else:
			rng=images
		for i in rng:
			tmp=numpy.dot(force[i].ravel(),force[i].ravel())
			tmp/=self.Atomcount
			ms+=tmp
		return numpy.sqrt(ms/self.numimages())



	def maxforce(self, images=None,force=None):
		"""return maximum force on single atom over images, use whole path if omitted
		@param images: if specified, only regard images in this list, otherwise use all images
		@param force: if specified, calculate rms of this forces array, otherwise use NEB forces
		"""
		if force==None:
			#TODO maxforce -> realforces
			force=self.nebforces
		maxF=0.0
		if images==None:
			rng=range(self.numimages())
		else:
			rng=images
		for i in rng:
			imforvec=numpy.reshape(force[i],(self.Atomcount,3))
			for j in range(self.Atomcount):
				tmpvec=imforvec[j]
				tmp=numpy.sqrt(numpy.dot(tmpvec,tmpvec))
				if tmp>maxF:
					maxF=tmp
		return maxF



	def writecenterdists(self, center=0, image=0, filename='neb.dst'):
		"""write distances of each atom from 'center' in 'image' into 'filename'
		@param center: center atom_bondcounts (default 0)
		@param image: image to calculate centerdists forcerms (default 0)
		@param filename: output filename
		"""
		outFile=open(filename,'w')
		for i in self.geos[image].centerdists(center):
			print("{0: 24.17E}".format(i),file=outFile)
		outFile.close



	def rmsds(self):
		"""return list of path RMSD w/r to image[0] of each atom"""
		rmsds=[]
		for i in range(self.Atomcount):
			tmp=0.0
			for j in range(1,self.numimages()):
				diff=self.geos[j].Geometry[i]-self.geos[0].Geometry[i]
				tmp+=numpy.dot(diff,diff)
			tmp/=(self.numimages()-1)
			rmsds.append(numpy.sqrt(tmp))
		return(rmsds)


	def writermsds(self, filename='rmsds.dat'):
		"""write path RMSD of each atom with. ref. to image[0] into output file
		@param filename: output filename
		"""
		ofile=open(filename,'w')
		for i in self.rmsds():
			print("{0: 24.17E}".format(i),file=ofile)
		ofile.close()




	def deltaforce(self, comparepath):
		"""return force difference vector F_self-F_comparepath
		comparepath path to calculate forces differences
		"""
		if (comparepath.Atomcount!=self.Atomcount or
			comparepath.numimages() != self.numimages()):
			raise("Paths must have same numbers of atoms and images")
		else:
			deltaf=[]
			for i in range(self.numimages()):
				deltaf.append(self.realforces[i]-comparepath.realforces[i])
			return deltaf



	def deltaenergy(self, comparepath):
		"""return total energy differences along self-comparepath
		@param comparepath: path to calculate energies difference
		"""
		if (comparepath.Atomcount!=self.Atomcount or
			comparepath.numimages() != self.numimages()):
			raise("Paths must have same numbers of atoms and images")
		else:
			deltaE=[]
			for i in range(self.numimages()):
				deltaE.append(self.energies[i]-comparepath.energies[i])
			return deltaE



	def splineResample(self, nnodes):
		""" resample geometries using natural cubic spline representation
		Attention! Splines will not be equally spaced along arc-length!
		L Populations will be deleted, Atom charges will be spline-interpolated
		@param nnodes: number of nodes to resample
		"""
		if self.verbosity>=constants.VBL_DEBUG2:
				print("Spline-resampling Reactionpath from nodecount")
		# construct list of evenly spaced parameter values
		positions=[float(i)/(nnodes-1.0) for i in range(nnodes)]
		# call splineListResample with the calulated positions
		self.splineListResample(positions)



	def splineListResample(self, paramlist):
		"""spline resample with nodes at specified parameter values
		@param paramlist: list of parameter values at which to place the new nodes"""
		if self.verbosity>=constants.VBL_DEBUG2:
				print("Spline-resampling Reactionpath from parameter list")
		# initialize spline representation, if necessary
		if not self.hasSplineRep:
			self._genSplineRep()
		# save important geometry information:
		geoshape=numpy.shape(self.geos[0].Geometry)
		AtomTypes=self.geos[0].AtomTypes
		Mode=self.geos[0].Mode
		Origin=self.geos[0].Origin
		AtomSubTypes=self.geos[0].AtomSubTypes
		AtomLayers=self.geos[0].AtomLayers
		LayerDict=self.geos[0].LayerDict
		# construct new geometries list
		newgeos=[]
		for position in paramlist:
			# interpolate data
			newcoordinates=numpy.reshape(self.splineRep["geo"].splint(position),geoshape)
			newlattice=numpy.reshape(self.splineRep["lat"].splint(position),(3,3))
			newcharges=list(self.splineRep["chr"].splint(position))
			# construct new image geometry object
			tempgeo=geometry.Geometry(Mode,self.Atomcount,AtomTypes,Origin,newlattice,newcoordinates,
				AtomLayers,LayerDict,newcharges,AtomSubTypes,None)
			# append new geometry to new geometries list
			newgeos.append(tempgeo)
		# now overwrite geometries arrays and reset energies and forces and kill old spline representation (for hygenic reasons)
		self.geos=newgeos
		#self.realforces=None
		#self.energies=[]
		self.splineRep=None
		self._rSplineRep=None



	def _genSplineRep(self):
		"""generate and store splien representation of path which interpolates coordinates, lattice and charges"""
		if self.verbosity>=constants.VBL_DEBUG2:
				print("generating spline representation of Reactionpath")
		nImages=self.numimages()
		parameter=numpy.zeros((nImages),dtype=float)
		vectors=numpy.zeros((nImages,self.Atomcount*3),dtype=float)
		lattices=numpy.zeros((nImages,9),dtype=float)
		charges=numpy.zeros((nImages,self.Atomcount),dtype=float)
		# construct vector spline nodes
		for i in range(nImages):
			parameter[i]=float(i)/float(nImages-1)
			vectors[i]=self.geos[i].Geometry.ravel()
			lattices[i]=self.geos[i].Lattice.ravel()
			charges[i]=numpy.array(self.geos[i].AtomCharges,dtype=float)
		# construct coordinate vectorSpline object
		geospline=spline.vectorSpline((parameter,vectors))
		# construct lattices vectorSpline object
		latticespline=spline.vectorSpline((parameter,lattices))
		# construct charges vectorSpline object
		chargespline=spline.vectorSpline((parameter,charges))
		self.splineRep={
			"geo":geospline,
			"lat":latticespline,
			"chr":chargespline
		}
	


	def getSplineArcLength(self):
		"""calculate total arc length of spline representation of Reactionpath
		@return: arc length of spline representation"""
		if not self.hasSplineRep:
			self._genSplineRep()
		return self.splineRep["geo"].arcLength
	
	splineArcLength=property(getSplineArcLength,doc="arc length of coordinates spline representation")




	def getSplineEquiParams(self, nparams):
		"""calculate a list spline representation parameter values that are arc-length equidistant in coordinates
		@param nparams: number of equidistant parameters to calculate
		@return: list of parameter values yielding arc-length equidistant coordinates from spline represetation"""
		arc=self.splineArcLength
		distance=arc/float(nparams-1)
		nsteps=int(1e4)
		arc=0.
		params=[0.]
		oldsample=self.splineRep["geo"].splint(0)
		for i in range(1,nsteps):
			newsample=self.splineRep["geo"].splint(float(i)*1e-4)
			diff=newsample-oldsample
			arc+=numpy.sqrt(numpy.dot(diff,diff))
			if arc >= distance:
				params.append(float(i)*1e-4)
				arc=0.
			oldsample=newsample
		params.append(1.)
		return params



	def equiDistSplineResample(self, nnodes):
		"""spline resample Reactionpath with equidistant nodes along spline
		@param nnodes: number of nodes in resampled path
		@return: distance between images"""
		#first get total arc length of spline representation
		params=self.getSplineEquiParams(nnodes)
		self.splineListResample(params)



	def _genRsplineRep(self):
		"""generate Renner-Subslpines representation of Path"""
		# first sanity-check: Renner Subspline Representation requires at least 5 images in path
		if self.numimages() < 5:
			raise(ValueError,"Attempt to generate Renner Subplines of path with less than 5 images.")
		# create copies of the flattened coordinate vectors at each image
		vectors=[]
		for i in range(self.numimages()):
			vectors.append(self.geos[i].Geometry.ravel())
		# initialize spline representation
		self._rSplineRep=spline.RennerSpline(vectors,closed=False)
		# finished



	def _updateRsplineRep(self):
		"""update Renner-Subspline Representation of self with current coordinates"""
		# create copies of the flattened coordinate vectors at each image
		vectors=[]
		for i in range(self.numimages()):
			vectors.append(self.geos[i].Geometry.ravel())
		# initialize spline representation
		self._rSplineRep.setNodes(vectors,closed=False)
		# finished



	def has_rSplineRep(self):
		"""check, if Path has Renner Subspline Representation
		@return: boolean, true if Renner Subplsine Representation has been set"""
		return self._rSplineRep!=None



	def rennerSubsplineResample(self,paramlist):
		"""resample path at lengths given in paramlist
		@param paramlist: Renner Subspline parameters (= lengths) at which to resample (default curve)
		"""
		# check if Renner Subpline Representation exists, bomb out if not, since parameter range depends on individual representation!
		if not self.has_rSplineRep():
			raise(ValueError,"Attempt to Renner-Subspline resample path but Renner Subspline representation not initialized")
		# save important geometry information:
		geoshape=numpy.shape(self.geos[0].Geometry)
		AtomTypes=self.geos[0].AtomTypes
		Mode=self.geos[0].Mode
		Origin=self.geos[0].Origin
		AtomSubTypes=self.geos[0].AtomSubTypes
		AtomLayers=self.geos[0].AtomLayers
		LayerDict=self.geos[0].LayerDict
		lattice=self.geos[0].Lattice
		# generate new geometries
		newgeos=[]
		for position in paramlist:
			# get interpolated coordinates
			newcoordinates=numpy.reshape(self._rSplineRep.splint(position),geoshape)
			# construct new image geometry object
			tempgeo=geometry.Geometry(Mode,self.Atomcount,AtomTypes,Origin,lattice,newcoordinates,
				AtomLayers,LayerDict,iAtomCharges=None,iAtomSubTypes=AtomSubTypes)
			# append new geometry to new geometries list
			newgeos.append(tempgeo)
		# now replace old path geometries and zero out related values
		self.geos=newgeos
		#self.energies=[]
		self.splineRep=None
		self._rSplineRep=None
		#self.realforces=None



	def movingAverageResample(self,windowsize):
		"""resample the path by averaging over windowsize images.
		New path has numimages-windowsize images
		@type windowsize: integer
		@param windowsize: number of images to average for each new image
		"""
		# first check if enough images are present for the requiested resampling
		if self.numimages() <= windowsize:
			raise ValueError("Not enough images present in reaction path to perform requested moving average")
		# build initial list of coordinate-vectors
		cVectors=["dummy"]
		for i in range(windowsize-1):
			cVectors.append(numpy.array(self.geos[i].Geometry.flat))
		# store some values for generation of new geometriy objects
		geoshape=self.geos[0].Geometry.shape
		geoclass=self.geos[0].__class__
		AtomTypes=self.geos[0].AtomTypes
		Mode=self.geos[0].Mode
		Origin=self.geos[0].Origin
		AtomSubTypes=self.geos[0].AtomSubTypes
		AtomLayers=self.geos[0].AtomLayers
		LayerDict=self.geos[0].LayerDict
		lattice=self.geos[0].Lattice
		# now iterate through new images
		newgeos=[]
		for i in range(windowsize-1,self.numimages()):
			del cVectors[0]
			cVectors.append(numpy.array(self.geos[i].Geometry.flat))
			newCoordinates=numpy.reshape(numpy.average(numpy.array(cVectors),0),geoshape)
			tempgeo=geoclass(Mode,self.Atomcount,AtomTypes,Origin,lattice,newCoordinates,
				iAtomLayers=AtomLayers,iLayerDict=LayerDict,iAtomCharges=None,iAtomSubTypes=AtomSubTypes)
			newgeos.append(tempgeo)
		# finished resampling, set new geometries etc
		self.geos=newgeos
		#self.energies=[]
		self.splineRep=None
		self._rSplineRep=None
		#self.realforces=None



	def hopCount(self, tolerance=0.1):
		"""
		count atomic position change events along the path.
		A position change is defines as any atomic coordinate changing by more than toleracne between two images
		@warning: This method is designed for use on vacancz position paths from ideal reference comparison and similar paths. It will absolute ly not work on MD trajectories.
		@param tolerance: position change detection threshold
		@return: total number of position change events along the path
		"""
		# Initialize
		lastgeo=self.geos[0].Geometry
		hopcounter=0
		# loop through all images except the first one
		for i in range(1,self.numimages()):
			# compare current Geometry to last image
			geoDifference=abs(self.geos[i].Geometry-lastgeo)
			# first create an per atom araray of the number of coordinates that changed by more than tolerance 
			temp1=numpy.add.reduce((geoDifference > 1),1)
			# now get hop count from counting all atoms which had mopre zan 0 coordinates change
			hopcounter+=numpy.add.reduce(temp1 > 0)
			# advance last image pointer
			lastgeo=self.geos[i].Geometry
		# finished, return
		return hopcounter

		

	def read_mdout(self, filename="md.out"):
		mdoutfile=utils.compressedopen(filename,"r")
		mdoutString="".join(list(mdoutfile))
		mdoutfile.close()
		self.parse_mdout(mdoutString)
		
	
	
	def parse_mdout(self, mdoutString):
		try:
			# read properties in atomic units from md.out
			mdoutLines = mdoutString.strip("\n").split("\n")
			# temperature
			temperatureLines = [i for i in mdoutLines if "temperature" in i.lower()]
			mdTemperatures = [float(i.split()[2]) for i in temperatureLines]
			# latticepressure
			pressureLines = [i for i in mdoutLines if "pressure" in i.lower()]
			mdPressure = [float(i.split()[1]) for i in pressureLines]
			# total energy
			etotLines = [i for i in mdoutLines if "total md energy" in i.lower()]
			mdEtot = [float(i.split()[3]) for i in etotLines]
			# ion potential energy
			epotLines = [i for i in mdoutLines if "potential energy" in i.lower()]
			mdEpot = [float(i.split()[2]) for i in epotLines]
			# ion kinetic energy
			ekinLines = [i for i in mdoutLines if "md kinetic energy" in i.lower()]
			mdEkin = [float(i.split()[3]) for i in ekinLines]		
			# md step
			stepLines = [i for i in mdoutLines if "md step" in i.lower()]
			mdStep = [int(i.split()[2]) for i in stepLines]
			# save properties
		except:
			print("Error while reading md.out file(s)")
			raise
		for image in xrange(0,self.numimages()):
			tempPos = mdStep.index(self.geos[image].timestep)
			if tempPos <= len(mdTemperatures):
				self.geos[image].iontemperature = mdTemperatures[tempPos]
				self.geos[image].latticepressure = mdPressure[tempPos]
				self.geos[image].totalenergy = mdEtot[tempPos]
				self.geos[image].ionpotentialenergy = mdEpot[tempPos]
				self.geos[image].ionkineticenergy = mdEkin[tempPos]
			else:
				print ("md.out file does not fit to the trajectory. It is too short!")
				raise
