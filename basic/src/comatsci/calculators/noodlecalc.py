##############################################################################
# noodlecalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
from .calculator import Calculator,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED#,CALCSTATUS_ERROR,CALCSTATUS_DISABLED

from .calcerror import CalcError
from .. import constants, utils
import ConfigParser

import tempfile
import os
import re
import shutil
import numpy



class noodlecalc(Calculator):
	"""Call NOODLE to calculate energies and forces"""


	defaults=dict(
		binary='noodle',
		skdir='SlKo',
		workdir='TEMP',
		chrdir='charges',
		rchr='true',
		paraminclude='params.ndl',
		infilename='dftb_in.hsd',
		oldSKnames='true'
		)


	def __init__(self, optionfname="pypath.ini", verbosity=1):
		"""construct NOODLE calculator
		@param optionfname: (default "pypath.ini") option file name
		@param verbosity: c.f. base class (default 1)		
		"""
		#@todo: replace option file name by passing a dictionary of configuration options
		Calculator.__init__(self, verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("initializing noodle calculator")
		# first parse config file and store into internal variables
		self.config = ConfigParser.SafeConfigParser(defaults=self.defaults)
		self.config.read(optionfname)
		if not self.config.has_section("NOODLE"):
			self.config.add_section("NOODLE")
		self.binary=self.config.get("NOODLE","binary")
		self.skdir=self.config.get("NOODLE","skdir")
		# if workdir directive=="TEMP" create a temporary directory as workdir
		self.workdir=self.config.get("NOODLE","workdir")
		if self.workdir=="TEMP":
			self.workdir=tempfile.mkdtemp(prefix="noodlecalc")
			self._rmworkdir=True
		else:
			self.workdir=os.path.abspath(self.workdir)
			if not os.path.exists(self.workdir):
				if self.verbosity>=constants.VBL_DEBUG1:
					print('noodle calculator: workdir "{0:s}" does not exist, creating it.'.format(self.workdir))()
				os.mkdir(self.workdir)
				self._rmworkdir=True
			else:
				self._rmworkdir=False
		self.chrdir=os.path.abspath(self.config.get("NOODLE","chrdir"))
		if not os.path.exists(self.chrdir):
			os.mkdir(self.chrdir)
		self.rchr=self.config.getboolean("NOODLE","rchr")
		self.paraminclude=self.config.get("NOODLE","paraminclude")
		self.infilename=self.config.get("NOODLE","infilename")
		self.remapatoms=None
		self.oldSKnames=self.config.getboolean("NOODLE","oldSKnames")



	def _worker(self):
		"""worker function to run noodle in a thread"""
		self.workreturncode=os.system(self.binary+" > noodle.out")
		return self.workreturncode



	def _writenoodleinput(self,Geo,charge,pchr=False):
		"""prepare the master NOODLE input file
		@param Geo: Geometry Object to perform calculation forget
		@param charge: total charge to pass to noodle
		@param pchr: should a point-charges file be used? (default False)		
"""
		#conversion from dftb to noodle maximum angular momentum spec
		maxang=["x","s","p","d","f"]
		ninput = open(self.infilename,"w")
		# first input the	user-provided parameters, then specify all our own data with override
		if self.paraminclude[-4:-1].lower==".xml":
			print("<<! {0:s}".format(self.paraminclude),file=ninput)
		else:
			print("<<+ {0:s}".format(self.paraminclude),file=ninput)
		#keep the geometry file separate
		print("""Geometry = GenFormat {\n <<< "input.gen" \n}""",file=ninput)
##		#this should hopefully give us forces and total energies at a single point
		#override initial charge reuse according to options	and existence of charge file
		if self.rchr:
			if os.path.exists("charges.bin"):
				print("*Hamiltonian = *DFTB {!ReadInitialCharges = Yes}",file=ninput)
		else:
			print("*Hamiltonian = *DFTB \{!ReadInitialCharges = No\}",file=ninput)
		#specify system charge
		print("*Hamiltonian = *DFTB {{!Charge = {0:f} }}".format(charge),file=ninput)
		#specify the Slater-Koster files and Max angular momenta
		sklist=[]
		mxalist=[]
		symlist,symdict=Geo.getatomsymlistdict() #@UnusedVariable
		for i in symlist:
			mxalist.append(Geo.PTE[i]+' = "'+maxang[Geo.LMAX[i]]+'"')
			for j in symlist:
				if self.oldSKnames:
					sklist.append(Geo.PTE[i]+"-"+Geo.PTE[j]+' = "./'
					+Geo.PTE[i].lower()+Geo.PTE[j].lower()+'"')
				else:
					sklist.append(Geo.PTE[i]+"-"+Geo.PTE[j]+' = "./'
					+Geo.PTE[i].capitalize()+"-"+Geo.PTE[j].capitalize()+'.skf"')
		newline="\n"
		print("*Hamiltonian = *DFTB {!SlaterKosterFiles = {",file=ninput)
		print(newline.join(sklist)+"}",file=ninput)
		print("!MaxAngularMomentum = {",file=ninput)
		print(newline.join(mxalist)+"}}",file=ninput)
		#override output options to our own needs
		print("""*Options = {
!AtomResolvedEnergies = No
!WriteResultsTag = Yes
!CalculateForces = Yes}""",file=ninput)
		#set pointcharges options, if specified
		if pchr:
			print("*Hamiltonian = *DFTB {*ElectricField ={ *PointCharges= {",file=ninput)
			print('<<< "pointcharges.xyzq"',file=ninput)
			print('}}}',file=ninput)
		ninput.close()



	def _prepare(self, steplabel, Geometry, charge):
		"""prepare NOODLE calculator run c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print("preparing noodle run")
		# if exitsts, copy old chages file
		chrfilename=steplabel+"-charges.bin"
		if os.path.exists(self.chrdir+"/"+chrfilename):
			shutil.copy(self.chrdir+"/"+chrfilename,self.rundir+"/charges.bin")
		# write the geometry files
		pchr=Geometry.layerbyname("PCHR")
		pchrenable= (pchr!=None)
		if pchrenable:
			Geometry.layersubgeometry(0).writegen("input.gen")
			self.remapatoms=(len(Geometry._layeratoms(0)),Geometry._layeratoms(0))
			Geometry.layersubgeometry(pchr).writexyzq("pointcharges.xyzq")
		else:
			Geometry.writegen("input.gen")
			#paranoia setting
			self.remapatoms=None
		self._writenoodleinput(Geometry,charge,pchrenable)
		# copy the SK files into the rundir
		symlist,symdict=Geometry.getatomsymlistdict() #@UnusedVariable
		for i in symlist:
			for j in symlist:
				if self.oldSKnames:
					shutil.copy(self.skdir+"/"+Geometry.PTE[i].lower()+Geometry.PTE[j].lower(),self.rundir)
				else:
					shutil.copy(self.skdir+"/"+Geometry.PTE[i].capitalize()+"-"+Geometry.PTE[j].capitalize()+".skf",self.rundir)
		# finally, copy the parameter file into the rundir
		shutil.copy(self.startdir+"/"+self.paraminclude,self.rundir)



	def _postrun(self, steplabel):
		"""Things to do after noodle run, i.e. save charges.bin, clean up c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print("noodle postrun cleanup and statistics")
		# first some statistics
		self.totalscf+=self.scfit
		self.totalruns+=1
		if self.verbosity>=constants.VBL_TALKY:
			print("{0:s}: SCC iterations: {1:3d}	 ----	 Total Energy: {2:12.6f} H".format(steplabel,self.scfit,self.etot))
		if os.path.exists(self.rundir+"/charges.bin"):
			chargefilename=steplabel+"-charges.bin"
			if not os.path.exists(self.chrdir):
				os.mkdir(self.chrdir)
			shutil.copy(self.rundir+"/charges.bin",self.chrdir+"/"+chargefilename)
		cleanuplist=os.listdir(".")
		for i in cleanuplist:
			os.unlink(i)
		os.chdir(self.startdir)
		os.rmdir(self.rundir)
		self.rundir=None



	def _readresults(self,atomcount):
		"""Read total energy and gradients from result files in current directory
		@param atomcount: number of atoms in system (ignored!)"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print("parsing noodle output")
		#overwrite atomcount, if we were working on a subgeometry
		if self.remapatoms!=None:
			realatomcount=atomcount
			atomcount=self.remapatoms[0]
		#first read the results.tag file into memory
		resultsfile=utils.compressedopen("results.tag","r")
		tagresults=utils.TaggedCollection(utils.ResultParser(resultsfile).entries)
		#get total energy
		self.etot=float(tagresults.getEntry("total_energy").value[0])
		#get forces
		if tagresults.getEntry("forces_calculated").value[0]:
			frcentry=tagresults.getEntry("forces")
			if frcentry.shape[1]!=atomcount:
				raise CalcError("number of forces returned by noodle does not match number of atoms")
			else:
				# funny reshaping construct to force a deep copy in order to get a contigous array
				temp=numpy.array(frcentry.value,dtype=float)
				self.gradients=numpy.reshape(temp,(atomcount,-1))
		else:
			raise CalcError("no noodle forces calculated")
		if tagresults.getEntry("scc").value[0]:
			self.scfit=int(tagresults.getEntry("n_scc_iters").value[0])
			if not tagresults.getEntry("scc_convergence").value[0] and self.verbosity>=constants.VBL_TALKY:
				print("Warning: noodle SCC not converged, forces may be wrong!")
		else:
			self.scfit=1
		#remap forces, in case we have been working on a subgeometry
		if self.remapatoms!=None:
			newgradients=numpy.zeros((realatomcount,3),dtype=float)
			for i in range(self.remapatoms[0]):
				newgradients[self.remapatoms[1][i]]=self.gradients[i]
			self.gradients=newgradients



	def shutdown(self):
		"""shut down the noodle calculator and delete workdir if it is a tempdir"""
		if self._rmworkdir:
			self.remove_workdir()
		Calculator.shutdown(self)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("noodle calculator shut down")
