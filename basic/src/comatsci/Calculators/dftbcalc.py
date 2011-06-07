##############################################################################
# dftbcalc.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from comatsci.Calculators.Calculator import *
from comatsci.Calculators.CalcError import *

try:
	from numpy.oldnumeric import *
except ImportError:
	from Numeric import *
	

class dftbcalc(Calculator):
	"""Class for dftb calculations and result retrieval"""

	#@deprecated: dftbcalc is tailored to Christof Koehler's version of the old dftb version. This version is no longer supported and not publicly available. Use noodlecalc with DFTB+ instead
	
	defaults=dict(
		binary="dftb",
		workdir="TEMP",
		chargesdir='charges',
		scc='t',
		econv='1e-6',
		cvonv='1e-3',
		rcharges='t',
		solver='1',		#0: ewevge, 1: dsygv, 2: dsyvd
		mixer='1',		#0: simple, 1: anderson, 2: broyden
		maxscc='35',
		amix='0.08',
		generations='4',
		keep_chrg='1',
		skdir='SK',
		tel='0.0',
		oldSKnames='t'
		)


	#Valence Electrons. semicore states are problematic!, this table ends after Ba!
	VALEL=[0,
	1,2,
	1,2,3,4,5,6,7,8,
	1,2,3,4,5,5,6,8,
	1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
	1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
	1,2]


	#lmax for dftb calculations.semicore states are problematic!, this table ends after Ba!
	LMAX=[0,
	1,1,
	2,2,2,2,2,2,2,2,
	3,3,3,3,3,3,3,3,
	4,4,4,4,4,4,4,4,4,4 ,4 ,4 ,4,4,4,4,4,4,
	5,5,5,5,5,5,5,5,5,5 ,5 ,5 ,5,5,5,5,5,5,
	6,6]


	def __init__(self, optionfname='pypath.ini', verbosity=1):
		"""initialite dftb calculator
		@param optionfname: (default) "pypath.ini" configuration file name
		@param verbosity: c.f. base class (default 1)
		"""
		#@todo: replace option file name by passing a dictionary of configuration options
		# Call the base class constructor
		Calculator.__init__(self,verbosity=verbosity)
		if self.verbosity>=constants.VBL_DEBUG1:
			print "initializing dftb calculator"
		# first parse config file and store into internal variables
		self.config = ConfigParser.SafeConfigParser(defaults=self.defaults)
		self.config.read(optionfname)
		self.binary=self.config.get("DFTB","binary")
		# if workdir directive=="TEMP" create a temporary directory as workdir
		self.workdir=self.config.get("DFTB","workdir")
		if self.workdir=="TEMP":
			self.workdir=tempfile.mkdtemp(prefix="dftbcalc")
			self._rmworkdir=True
		else:
			self.workdir=os.path.abspath(self.workdir)
			if not os.path.exists(self.workdir):
				if self.verbosity>=constants.VBL_DEBUG1:
					print 'dftb calculator: workdir "%s" does not exist, creating it.'
				os.mkdir(self.workdir)
				self._rmworkdir=True
			else:
				self._rmworkdir=False
		self.chargesdir=os.path.abspath(self.config.get("DFTB","chargesdir"))
		self.skdir=os.path.abspath(self.config.get("DFTB","skdir"))
		self.rcharges=self.config.get("DFTB","rcharges")
		self.scc=self.config.get("DFTB","scc")
		self.solver=self.config.get("DFTB","solver")
		self.mixer=self.config.get("DFTB","mixer")
		self.maxscc=self.config.get("DFTB","maxscc")
		self.econv=self.config.get("DFTB","econv")
		self.cconv=self.config.get("DFTB","cconv")
		self.generations=self.config.get("DFTB","generations")
		self.tel=self.config.get("DFTB","tel")
		self.amix=self.config.get("DFTB","amix")
		self.scc=self.config.get("DFTB","scc")
		self.keep_chrg=self.config.get("DFTB","keep_chrg")
		self.oldSKnames=self.config.getboolean("DFTB","oldSKnames")



	def _writedftbinput(self,electrons,atomcount,symlist,symdict):
		"""write dftb input file called dftb.in in current directory
		@param electrons: number of valence electrons in system
		@param atomcount: number of atoms in system
		@param symlist: list of element symols in system
		@param smydict: dictionary of (atom type number:	element symbol)-s"""
		inp=open("dftb.in","w")
		#first dftb input line
		print >> inp, "3 0.001 "+self.scc+" "+self.econv+" ",
		print >> inp, self.cconv+" ",
		if not os.path.exists("CHR.DAT"):
			print >> inp, "F F"
		else:
			print >> inp, self.rcharges+" F"
		#second dftb input line
		print >> inp, "'input.gen'"
		#number of electrons
		print >> inp, electrons
		#Solver, Mixer, ... line
		print >> inp, self.solver+" "+ self.mixer+" "+ self.maxscc+" "+ self.amix+" ",
		print >> inp, self.generations+" "+self.keep_chrg+" 1"
		#Number of movable atoms
		print >> inp, atomcount
		#output geometry name (irrelevant)
		print >> inp, "'output.gen'"
		#now lmax by type index
		for i in symlist:
			print >> inp, str(self.LMAX[i])+" ",
		print >> inp
		#now the SK file list
		if self.oldSKnames:
			for i in symlist:
				for j in symlist:
					print >> inp,"'"+symdict[i].lower()+symdict[j].lower()+"'"
		else:
			for i in symlist:
				for j in symlist:
					print >> inp,"'"+symdict[i].capitalize()+"-"+symdict[j].capitalize()+".skf'"
		#final input line, only t_el is of interest...
		print >> inp, "40 "+self.tel+" "+self.tel+" "+"0.2 0"
		inp.close()


	def _worker(self):
		"""worker function to run dftb in a thread"""
		self.workreturncode=os.system(self.binary+" < dftb.in > dftb.out")
		return self.workreturncode


	def _postrun(self, steplabel):
		"""Things to do after dftb run, i.e. save CHR.DAT, clean up c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print "dftb postrun cleanup and statistics"
		# first some statistics
		self.totalscf+=self.scfit
		self.totalruns+=1
		if self.verbosity>=constants.VBL_TALKY:
			print "%s: SCC iterations: %3d   ----   Total Energy: %12.6fH" % (steplabel,self.scfit,self.etot)
		chargefilename=steplabel+"-CHR.DAT"
		if not os.path.exists(self.chargesdir):
				os.mkdir(self.chargesdir)
		shutil.copy(self.rundir+"/CHR.DAT",self.chargesdir+"/"+chargefilename)
		cleanuplist=os.listdir(".")
		for i in cleanuplist:
			os.unlink(i)
		os.chdir(self.startdir)
		os.rmdir(self.rundir)
		self.rundir=None



	def _readresults(self,atomcount):
		"""Read total energy and gradients from result files in current directory
		@param atomcount: number of atoms in system"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print "parsing dftb output"
		if self.status()!=CALCSTATUS_FINISHED:
			raise CalcError('Try to get results form unfinished calculation')
		# first check for diagonalization errors
		if os.path.exists("dftb.out") or os.path.exists("dftb.out.gz") or os.path.exists("dftb.out.bz2"):
			ofile=utils.compressedopen("dftb.out","r")
			lines=list(ofile)
			dummy=lines[-1].split()
			if dummy[0].lower()=='ewevge:':
				raise CalcError("Diagon Error")
			else:
				dummy=lines[-3].split()
				self.workercputime+=float(dummy[2])
		else:
			raise CalcError("Results unreadable")
		# Check if parseable output exists
		if (not os.path.exists("ENERGY.TMP") or os.path.exists("ENERGY.TMP.gz") or os.path.exists("ENERGY.TMP.bz2")) or (not os.path.exists("FRC.DAT") or os.path.exists("FRC.DAT.gz") or os.path.exists("FRC.DAT.bz2")):
				raise CalcError("Results unreadable")
		# now read the energy
		enfile=utils.compressedopen("ENERGY.TMP")
		for i in range(4):
			line=enfile.readline()
		dummy=line.split()
		self.scfit=int(dummy[2])
		self.etot=float(dummy[4])
		enfile.close
		if self.scfit>=int(self.maxscc) and self.verbosity>=constants.VBL_QUIET:
			print " *** WARNING: SCC probably not converged"
			print " maxscc : %s scc iterations : %d" % (self.maxscc,self.scfit)
		gradfile=utils.compressedopen("FRC.DAT")
		gradsbuf=[]
		for i in range(5):
			line=gradfile.readline()
		for i in range(atomcount):
			line=gradfile.readline()
			dummy=line.split()
			# This will break, if forces are not ordered in FRC.DAT!
			if dummy[0][0].isdigit():
				gradsbuf.append([ float(s) for s in dummy[1:4] ])
			else:
				for j in range(i,atomcount):
					gradsbuf.append([ 0.0, 0.0, 0.0 ])
				break
		gradfile.close()
		self.gradients=array(gradsbuf)


	def _prepare(self, steplabel, Geometry, charge):
		"""prepare dftb calculator run c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print "preparing dftb run"
		# if exitsts, copy old charge file
		chargefilename=steplabel+"-CHR.DAT"
		if os.path.exists(self.chargesdir+"/"+chargefilename):
			shutil.copy(self.chargesdir+"/"+chargefilename,self.rundir+"/CHR.DAT")
		# write the geometry file
		Geometry.writegen("input.gen")
		#now write dftb.input file
		symlist,symdict=Geometry.getatomsymlistdict()
		self._writedftbinput(Geometry.velcount(self.VALEL)-charge,
			Geometry.Atomcount,symlist,Geometry.PTE)
		# finally, copy the SK files into the rundir
		for i in symlist:
			for j in symlist:
				if self.oldSKnames:
					shutil.copy(self.skdir+"/"+Geometry.PTE[i].lower()+Geometry.PTE[j].lower(),self.rundir)
				else:
					shutil.copy(self.skdir+"/"+Geometry.PTE[i].capitalize()+"-"+Geometry.PTE[j].capitalize()+".skf",self.rundir)


	def runfg(self, Geometry, steplabel,charge=0):
		"""run calculation in foreground, try ewevge if dsygv* failed c.f. base class"""
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



	def shutdown(self):
		"""shut down the dftb calculator and deletw workdir if it is a tempdir"""
		if self._rmworkdir:
			self.remove_workdir()
		Calculator.shutdown(self)
		if self.verbosity>=constants.VBL_DEBUG1:
			print "dftb calculator shut down"

