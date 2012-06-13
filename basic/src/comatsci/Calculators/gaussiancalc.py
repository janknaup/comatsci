##############################################################################
# gaussiancalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function
from comatsci.Calculators.Calculator import Calculator#,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED,CALCSTATUS_ERROR,CALCSTATUS_DISABLED

from comatsci.Calculators.CalcError import CalcError
import comatsci.constants as constants
import comatsci.utils as utils
import ConfigParser
import tempfile
import os
#import sys
import shutil
import numpy
#import time
import re	

class gaussiancalc(Calculator):
	"""Class for Gaussian (03) calculations and result retrieval"""
	#TODO: add point charges support to gaussiancalc

	defaults=dict(
		binary='g03',
		workdir='TEMP',
		chkdir='g03chk',
		rchk='false',
		hamiltonian='',
		link0lines='',
		routeopts='',
		spinmul="1"
		)


	def __init__(self, optionfname="pypath.ini", verbosity=1):
		"""construct Gaussian (03) calculator
		@param optionfname: (default "pypath.ini") option file name
		@param verbosity: c.f. base class (default 1)
		"""
		#@todo: replace option file name by passing a dictionary of configuration options
		Calculator.__init__(self, verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("initializing gaussien(03) calculator")
		# first parse config file and store into internal variables
		self.config = ConfigParser.SafeConfigParser(defaults=self.defaults)
		self.config.read(optionfname)
		if not self.config.has_section("Gaussian"):
			self.config.add_section("Gaussian")
		self.binary=self.config.get("Gaussian","binary")
		# if workdir directive=="TEMP" create a temporary directory as workdir
		self.workdir=self.config.get("Gaussian","workdir")
		if self.workdir=="TEMP":
			self.workdir=tempfile.mkdtemp(prefix="gaussiancalc")
			self._rmworkdir=True
		else:
			self.workdir=os.path.abspath(self.workdir)
			if not os.path.exists(self.workdir):
				if self.verbosity>=constants.VBL_DEBUG1:
					print('gaussian calculator: workdir "{0:s}" does not exist, creating it.')
				os.mkdir(self.workdir)
				self._rmworkdir=True
			else:
				self._rmworkdir=False
		self.chkdir=os.path.abspath(self.config.get("Gaussian","chkdir"))
		if not os.path.exists(self.chkdir):
			os.mkdir(self.chkdir)
		self.rchk=self.config.getboolean("Gaussian","rchk")
		self.hamiltonian=self.config.get("Gaussian","hamiltonian")
		if self.hamiltonian=='':
			raise CalcError("No hamiltonian specified! abort")
		self.link0lines=self.config.get("Gaussian","link0lines")
		self.routeopts=self.config.get("Gaussian","routeopts")
		self.spinmul=self.config.getint("Gaussian","spinmul")
		if self.verbosity>=constants.VBL_DEBUG1:
			print("gaussian calculator initialized")


	def _writegaussianinput(self, steplabel, Geometry, charge=0.0):
		"""prepare the Gaussian input file
		@param steplabel: name of current calculation
		@param charge: system total charge"""
		ginput=open('input.com','w')
		#Link 0 section
		print(" %chk=path.chk",file=ginput)
		if self.link0lines!='':
			print(self.link0lines,file=ginput)
		#Route section
		routeline=" #"
		routeline+=self.hamiltonian+" FChk=ForceCart"
		routeline+=" force "
		routeline+=" "+self.routeopts
		print(routeline,file=ginput)
		print(" ",file=ginput)
		#Title section
		print(steplabel+"\n",file=ginput)
		#Molecule spec
		print(" {0:d} {1:d}".format(charge,self.spinmul),file=ginput)
		print(Geometry.gaussianstring(),file=ginput)
		#!trailing empty line!
		print("\n",file=ginput)



	def _worker(self):
		"""worker function to run Gaussian"""
		self.workreturncode=os.system(self.binary+" < input.com > output.g03")
		return self.workreturncode


	def _postrun(self, steplabel):
			"""Things to do after Gaussian run, i.e. save .com file, clean up
			@param steplabel: name of current calculation"""
			if self.verbosity>=constants.VBL_DEBUG2:
				print("gaussian postrun cleanup and statistics")
			# first some statistics
			self.totalscf+=self.scfit
			self.totalruns+=1
			if self.verbosity>=constants.VBL_TALKY:
				print("{0:s}: SCF iterations: {1:3d}   ----   Total Energy: {2:12.6f} H".format(steplabel,self.scfit,self.etot))
			chkfilename=steplabel+".chk"
			if not os.path.exists(self.chkdir):
				os.mkdir(self.chkdir)
			shutil.copy(self.rundir+"/path.chk",self.chkdir+"/"+chkfilename)
			cleanuplist=os.listdir(".")
			for i in cleanuplist:
				os.unlink(i)
			os.chdir(self.startdir)
			os.rmdir(self.rundir)
			self.rundir=None



	def _readresults(self, atomcount):
		"""prese results from Gaussian(03) calculation
		@param atomcount: Number of atoms in calculation"""
		if self.verbosity>=constants.VBL_DEBUG1:
			print("parsing gaussian output")
		# read fchk file into Memory
		fchkfile=utils.compressedopen("Test.FChk")
		fchklines=list(fchkfile)
		fchkfile.close()
		# compile regular expressions for Total Energy and forces
		energyre=re.compile("^Total Energy")
		forcere=re.compile("^Cartesian Gradient")
		#iterate through lines
		for i in range(len(fchklines)):
			nrgparsed=False
			frcparsed=False
			if not nrgparsed and energyre.search(fchklines[i]):
				dummy=fchklines[i].strip().split()
				self.etot=float(dummy[-1])
				nrgparsed=True
			if not frcparsed and forcere.search(fchklines[i]):
				dummy=fchklines[i].strip().split()
				parsedcomponents=0
				numcomponents=int(dummy[-1])
				if numcomponents!=atomcount*3:
					raise CalcError("Returned number of force components does not match number of atoms.")
				currentline=i+1
				tempgrad=[]
				while parsedcomponents<numcomponents:
					dummy=fchklines[currentline].split()
					for j in dummy:
						tempgrad.append(float(j))
						parsedcomponents+=1
					currentline+=1
				self.gradients=numpy.reshape(numpy.array(tempgrad,dtype=float),(-1,3))
				frcparsed=True
			if frcparsed and nrgparsed:
				break



	def _prepare(self, steplabel, Geometry, charge):
		"""prepare Gaussian calculator run c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG1:
			print("preparing gaussian run")
		# if exitsts and we want to read it, copy old DM file
		chkfilename=steplabel+".chk"
		if os.path.exists(self.chkdir+"/"+chkfilename) and self.rchk:
			shutil.copy(self.chkdir+"/"+chkfilename,self.rundir+"/path.chk")
		# write the input file
		self._writegaussianinput(steplabel,Geometry, charge)




	def shutdown(self):
		"""shut down the siesta calculator and deletw workdir if it is a tempdir"""
		if self._rmworkdir:
			self.remove_workdir()
		Calculator.shutdown(self)
		if self.verbosity>=constants.VBL_DEBUG1:
			print("gaussian calculator shut down")
