##############################################################################
# siestacalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
from comatsci.Calculators.Calculator import Calculator,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED#,CALCSTATUS_ERROR,CALCSTATUS_DISABLED

from comatsci.Calculators.CalcError import CalcError
import comatsci.constants as constants
import comatsci.utils as utils
import ConfigParser
import tempfile
import os
#import sys
import shutil
import numpy
import re	


class siestacalc(Calculator):
	"""Class for SIETSA calculations and result retrieval"""

	defaults=dict(
		binary='siesta',
		ppdir='pseudo',
		workdir='TEMP',
		dmdir='DMAPS',
		paraminclude='params.fdf'
		)


	def __init__(self, optionfname="pypath.ini", verbosity=1):
		"""construct SIETSA calculator
		@param optionfname: (default "pypath.ini") option file name
		@param verbosity: c.f. base class (default 1)
		"""
		#@todo: replace option file name by passing a dictionary of configuration options
		Calculator.__init__(self, verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("initializing SIESTA calculator")
		# first parse config file and store into internal variables
		self.config = ConfigParser.SafeConfigParser(defaults=self.defaults)
		self.config.read(optionfname)
		if not self.config.has_section("SIESTA"):
			self.config.add_section("SIESTA")
		self.binary=self.config.get("SIESTA","binary")
		self.ppdir=self.config.get("SIESTA","ppdir")
		# if workdir directive=="TEMP" create a temporary directory as workdir
		self.workdir=self.config.get("SIESTA","workdir")
		if self.workdir=="TEMP":
			self.workdir=tempfile.mkdtemp(prefix="siestacalc")
			self._rmworkdir=True
		else:
			self.workdir=os.path.abspath(self.workdir)
			if not os.path.exists(self.workdir):
				if self.verbosity>=constants.VBL_DEBUG1:
					print('siesta calculator: workdir "{0:s}" does not exist, creating it.'.format(self.workdir))
				os.mkdir(self.workdir)
				self._rmworkdir=True
			else:
				self._rmworkdir=False
		self.dmdir=os.path.abspath(self.config.get("SIESTA","dmdir"))
		if not os.path.exists(self.dmdir):
			os.mkdir(self.dmdir)
		self.rdms=self.config.getboolean("SIESTA","rdms")
		self.paraminclude=self.config.get("SIESTA","paraminclude")
		if self.verbosity>=constants.VBL_DEBUG1:
			print("siesta calculator initialized")


	def _writesiestainput(self, steplabel, charge=0.0):
		"""prepare the master SIESTA input file
		@param steplabel: name of current calculation
		@param charge: system total charge"""
		sinput = open("input.fdf","w")
		print("SystemLabel path",file=sinput)
		print("SystemName comatsci calculation for {0:s}".format(steplabel),file=sinput)
		if self.rdms:
			if os.path.exists("path.DM"):
				print("DM.useSaveDM .true.",file=sinput)
			else:
				print("DM.useSaveDM .false.",file=sinput)
		else:
			print("DM.useSavedDM .false.",file=sinput)
		print("MD.TypeOfRun CG\nMD.NumCGSteps 0",file=sinput)
		print("NetCharge {0:d}".float(float(charge)),file=sinput)
		print("%include {0:s}\n%include geometry.fdf".format(self.paraminclude),file=sinput)
		sinput.close()



	def _worker(self):
		"""worker function to run SIESTA"""
		self.workreturncode=os.system(self.binary+" < input.fdf > output.fdf")
		return self.workreturncode


	def _postrun(self, steplabel):
			"""Things to do after SIESTA run, i.e. save DM file, clean up
			@param steplabel: name of current calculation"""
			if self.verbosity>=constants.VBL_DEBUG2:
				print("SIESTA postrun statistics and cleanup")
			# first some statistics
			self.totalscf+=self.scfit
			self.totalruns+=1
			if self.verbosity>=constants.VBL_TALKY:
				print("{0:s}: SCF iterations: {1:3d}   ----   Total Energy: {2:12.6f} H".format(steplabel,self.scfit,self.etot))
			dmfilename=steplabel+".DM"
			if not os.path.exists(self.dmdir):
				os.mkdir(self.dmdir)
			shutil.copy(self.rundir+"/path.DM",self.dmdir+"/"+dmfilename)
			cleanuplist=os.listdir(".")
			for i in cleanuplist:
				os.unlink(i)
			os.chdir(self.startdir)
			os.rmdir(self.rundir)
			self.rundir=None



	def _readresults(self,atomcount):
			"""Read total energy and gradients from result files in current directory
			@param atomcount: number of atoms in this calculation"""
			if self.verbosity>=constants.VBL_DEBUG2:
				print("parsing SIESTA results")
			if self.status()!=CALCSTATUS_FINISHED:
				raise CalcError('Try to get results form unfinished calculation')
			# Check if parseable output exists
			if (not os.path.exists("output.fdf") or os.path.exists("output.fdf.gz") or os.path.exists("output.fdf.bz2")) or (not os.path.exists("path.FA") or os.path.exists("path.FA.gz") or os.path.exists("path.FA.bz2")):
					raise CalcError("Results unreadable")
			# now read the energy
			tere=re.compile("^siesta:.*Total =")
			itre=re.compile("timer:  IterSCF")
			enfile=utils.compressedopen("output.fdf")
			for line in enfile:
				if tere.search(line):
					dummy=line.split()
					self.etot=float(dummy[3])*constants.EVOLT
				if itre.search(line):
					dummy=line.split()
					self.scfit=int(dummy[2])-1
			enfile.close
			gradfile=utils.compressedopen("path.FA")
			gradsbuf=[]
			atomcount=int(gradfile.readline())
			for i in range(atomcount):
				line=gradfile.readline()
				dummy=line.split()
				# This will break, if forces are not ordered in FRC.DAT!
				if dummy[0][0].isdigit():
					gradsbuf.append([ (float(s)*constants.EVPERANG) for s in dummy[1:4] ])
			gradfile.close()
			self.gradients=numpy.array(gradsbuf)


	def _prepare(self, steplabel, Geometry, charge):
		"""prepare SIESTA calculator run c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG1:
			print("preparing SIESTA run")
		# if exitsts, copy old DM file
		dmfilename=steplabel+".DM"
		if os.path.exists(self.dmdir+"/"+dmfilename):
			shutil.copy(self.dmdir+"/"+dmfilename,self.rundir+"/path.DM")
		# write the geometry file
		Geometry.writefdf("geometry.fdf")
		self._writesiestainput(steplabel,charge)
		# copy the parameter file into the rundir
		shutil.copy(self.startdir+"/"+self.paraminclude,self.rundir)
		# finally, copy the PP files into the rundir
		symlist,symdict = Geometry.getatomsymlistdict()
		for i in symlist:
			if os.path.exists(self.ppdir+"/"+Geometry.PTE[i]+".psf"):
				shutil.copy(self.ppdir+"/"+Geometry.PTE[i]+".psf",self.rundir)
			elif os.path.exists(self.ppdir+"/"+Geometry.PTE[i]+".psf.gz"):
				shutil.copy(self.ppdir+"/"+Geometry.PTE[i]+".psf.gz",self.rundir)
				os.system("gzip -d "+self.rundir+"/"+Geometry.PTE[i]+".psf.gz")
			else:
				raise CalcError("SIESTA PP file not found")



	def shutdown(self):
		"""shut down the siesta calculator and deletw workdir if it is a tempdir"""
		if self._rmworkdir:
			self.remove_workdir()
		Calculator.shutdown(self)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("siesta calculator shut down")
