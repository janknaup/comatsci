##############################################################################
# Calculator.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

try:
	from numpy.oldnumeric import *
except ImportError:
	from Numeric import *

import os
import sys
import shutil
import ConfigParser
import time
import tempfile
import re
from comatsci import Spline
from comatsci import constants
from comatsci import utils
from comatsci.Calculators import CalcError

# Calculator Status Flags
CALCSTATUS_READY=0
CALCSTATUS_RUNNING=1
CALCSTATUS_FINISHED=2
CALCSTATUS_ERROR=-1
CALCSTATUS_DISABLED=-2

# Some useful constants
# This converts eV to H
EVOLT=constants.EVOLT
#HARTREE=27.211384523 H to eV conversion
HARTREE=constants.HARTREE
BOHR=constants.BOHR
ANGSTROM=constants.ANGSTROM
# This converts eV/Ang to H/Bohr
EVPERANG=constants.EVPERANG


class Calculator:

	"""Base Class to do total Energy and forces calculations"""
	def __init__(self, verbosity=None):
		"""This constructor should be called within the constructor of any derived class!
		@param verbosity: level of output verbosity if unspecified, choose VBL_NORMAL (default None)
		"""
		# initialize some work variables
		self.startdir=os.getcwd()
		self.etot=None
		self.gradients=None
		self.scfit=0
		if verbosity==None:
			self.verbosity=constants.VBL_NORMAL
		else:
			self.verbosity=verbosity
		# some statistics for profiling
		self.totalscf=0
		self.workercputime=0.
		self.workerwalltime=0.
		self.totalruns=0
		self.calculatorcputime=0.
		self.calculatorwalltime=0.
		# book keeping
#		self.WORKTHREAD=threading.Thread(target=self._worker,name="WORKER")
		self.rundir=None
		self.workreturncode=None
		self._status=CALCSTATUS_READY



	def remove_workdir(self):
		"""this method removes the workdir, in order to clean up temporary directories"""
		# delete the working directory if it still exists. CAVEAT!
		if os.path.exists(self.workdir):
			os.chdir(self.workdir+"/..")
			shutil.rmtree(self.workdir)



	def status(self):
		"""return status of calculator. check the 
		calculator thread's vitality and exit code"""
		if self._status != CALCSTATUS_RUNNING:
			return self._status
#		elif self.WORKTHREAD.isAlive():
#			return CALCSTATUS_RUNNING
		elif self.workreturncode==0:
			self._status=CALCSTATUS_FINISHED
			return CALCSTATUS_FINISHED
		else:
			self._status=CALCSTATUS_ERROR
			return CALCSTATUS_ERROR



	def getforces(self):
		"""return an N*3 array of the forces"""
		if self.status()==CALCSTATUS_FINISHED:
			return self.gradients
		else:
			raise CalcError("Calculator not finished")



	def getenergy(self):
		"""return the total energy"""
		if self.status()==CALCSTATUS_FINISHED:
			return self.etot
		else:
			raise CalcError("Calculator not finished")



	def _worker(self):
		"""DUMMY! - This is the worker function of the WORKTHREAD
		it should contain the system call (or whatever) to run the actual binary"""
		pass



	def _prepare(self, steplabel, Geometry, Charge):
		"""DUMMY! - This method is called directly after entering the rundir and
		before starting the WORKTHREAD. It should create all necessary input files
		best implement this to be independent of the current path...
		@param steplabel: name to use for current calculation
		@param Geometry: Geometry object to calculate for
		@param Charge: total charge to pass to calcultion
		"""
		pass


	def _readresults(self):
		"""DUMMY! - This method is meant to read the results of the calculation and
		store energy and forces in memory"""


	def _postrun(self,steplabel):
		"""DUMMY! - This method is meant to clean up after the calculation is done,
		to store store any restart files and to delete the rundir when finished
		@param steplabel: name of current calculation
		"""


	def start(self, Geometry, steplabel,charge=0):
		"""Start WORKTHREAD
		@param Geometry: Geometry object to calculate forcerms
		@param steplabel: name of current calculation
		@param charge: total charge to calculate with (default 0)
		"""
		if self.status() != CALCSTATUS_READY:
			raise CalcError('Calculator not ready')
		self.rundir=self.workdir+"/"+steplabel
		if not os.path.exists(self.rundir):
			os.mkdir(self.rundir)
		os.chdir(self.rundir)
		self._prepare(steplabel, Geometry, charge)
		#now spawn dftb thread
		#this will probably only work under *NIX, but who cares?
#		self.WORKTHREAD.start()
		if self.verbosity>=constants.VBL_TALKY:
			print "started calculation for step: %s." %(steplabel)
		self._status=CALCSTATUS_RUNNING



	def runfg(self, Geometry, steplabel,charge=0):
		"""run calculation in foreground
		@param Geometry: Geometry object to calculate for
		@param steplabel: name of current calculation
		@param charge: total charge to calculate with (default 0)
		"""
		error=None
		startcpu=time.clock()
		startwall=time.time()
		if self.status() != CALCSTATUS_READY:
			raise CalcError('Calculator not ready')
		self.rundir=self.workdir+"/"+steplabel
		if not os.path.exists(self.rundir):
			os.mkdir(self.rundir)
		os.chdir(self.rundir)
		self._prepare(steplabel, Geometry, charge)
		if self.verbosity>=constants.VBL_TALKY:
			print "running calculation for step: %s." %(steplabel)
		self._status=CALCSTATUS_RUNNING
		self._worker()
		if self.verbosity>=constants.VBL_TALKY:
			print "calculation done for %s" % (steplabel)
		# read the results
		self._readresults(Geometry.Atomcount)
		# clean up
		self._postrun(steplabel)
		self.calculatorcputime+=(time.clock()-startcpu)
		self.calculatorwalltime+=(time.time()-startwall)



	def calcInput(self, Geometry, steplabel,charge=0):
		"""Just prepare the calculation and leave job directory as it is
		@param Geometry: Geometry object to calculate for
		@param steplabel: name of current calculation
		@param charge: total charge to calculate with (default 0)
		"""
		#store path prior to input generation
		priorpath=os.path.realpath(".")
		#avoid interference with running calculations or losing unprocessed results
		if self.status() != CALCSTATUS_READY:
			raise CalcError('Calculator not ready')
		#prepare input
		self.rundir=self.workdir+"/"+steplabel
		if not os.path.exists(self.rundir):
			os.mkdir(self.rundir)
		os.chdir(self.rundir)
		if self.verbosity>=constants.VBL_TALKY:
			print "preparing run for step: %s." %(steplabel)
		self._prepare(steplabel, Geometry, charge)
		#return to path prior to input generation
		os.chdir(priorpath)



	def parseForeignOutput(self, atomCount, dataDir):
		"""Parse output from calculator program to use results from external calculations insite comatsci<br>
		@param atomCount: Number of atoms in Geometry
		@param dataDir: directory from which to reead the calculator output
		"""
		priorpath=os.path.realpath(".")
		os.chdir(dataDir)
		self._status=CALCSTATUS_FINISHED
		self._readresults(atomCount)
		os.chdir(priorpath)



	def finreready(self):
		"""re-ready the calculator, after finished calculation
		obliterating all results in memory"""
		if self._status==CALCSTATUS_FINISHED:
			self.etot=None
			self.gradients=None
#			self.WORKTHREAD=threading.Thread(target=self._worker,name="WORKER")
			self._status=CALCSTATUS_READY
			self.scfit=0
			self.workreturncode=None
		else:
			raise CalcError("Calculator not finished")



	def errorreready(self):
		"""re-ready the calculator, after an error occured
		obliterating all results in memory"""
		if self._status==CALCSTATUS_ERROR:
			self.etot=None
			self.gradients=None
#			self.WORKTHREAD=threading.Thread(target=self._worker,name="WORKER")
			self._status=CALCSTATUS_READY
			self.scfit=0
			self.workreturncode=None
			# clean up directory of broken run...
			if (self.rundir!=None) and os.path.exists(self.rundir):
				os.chdir(self.rundir)
				cleanuplist=os.listdir(".")
				for i in cleanuplist:
					os.unlink(i)
				os.chdir(self.startdir)
				os.rmdir(self.rundir)
				self.rundir=None
		else:
			raise CalcError("Error Cleanup Routine Called, but no error state")



	def getiterations(self):
		"""return number of iterations until convergence"""
		return self.scfit
	iterations=property(getiterations,doc="the number of iteratons performed in the last run")
	
	
	
	def run(self, options):
		"""run the calculator in the foreground, 
		unpacking all arguments from the options dictionary.
		rereadies the calculator in if it's status is FINISHED, to hide status
		tracking from schedulers.
		For use with generic schedulers.
		@param : options Dictionary containing (at least) a string "steplabel", a geometry 
		object "Geometry" and a float "Charge" corresponding to the arguments
		of runfg()
		"""
		if self._status==CALCSTATUS_FINISHED:
			self.finreready()
		self.runfg(options["Geometry"],options["steplabel"],options["Charge"])
		
		
	
	def getresults(self):
		"""return the results of the last calculation as a dictionary.
		for use with generic schedulers.
		"""
		#deepcopy the results, to make sure no implementation of
		#finreready kills out return values
		resultsdict={
			"Forces":copy.deepcopy(self.getforces()),
			"Energy":copy.deepcopy(self.getenergy())
			}
		return resultsdict



	def shutdown(self):
		"""shut down the claculator"""
		self._status=CALCSTATUS_DISABLED
