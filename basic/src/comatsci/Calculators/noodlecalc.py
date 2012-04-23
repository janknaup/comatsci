##############################################################################
# noodlecalc.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from comatsci.Calculators.Calculator import Calculator,CALCSTATUS_READY,CALCSTATUS_RUNNING,CALCSTATUS_FINISHED#,CALCSTATUS_ERROR,CALCSTATUS_DISABLED

from comatsci.Calculators.CalcError import CalcError
import comatsci.constants as constants
import comatsci.utils as utils
import ConfigParser

import tempfile
import os
import re
import shutil
import numpy


################################################################################
#  Helper classes to robustly parse noodle output
#  written by Balint Aradi <Aradi@bccms.uni-bremen.de>, copy & pasted by Jan Knaup
################################################################################

############################################################################
# Conversion functors
############################################################################

class ConversionError(Exception):
	pass

class InvalidEntry(Exception):
	pass

class Converter(object):
  """Base class for string to value converters"""

  def __init__(self, nolist=False):
    """nolist -- if True, not a list, but only a single value is returned after
                 conversion.
    """
    self.nolist = nolist

  def __call__(self, strValue):
    """strValue -- string represenation of the values to convert"""
    values = strValue.split()
    if self.nolist and len(values) > 1:
      raise ConversionError, "Too many values"
    result = self.convert(values)
    if self.nolist:
      return result[0]
    else:
      return result

  def convert(self, values):
    """Conversion function
    values -- list of strings representation of values to convert
    """
    return values



class FloatConverter(Converter):
  """Converts string to float"""

  def convert(self, values):
    ll = []
    for val in values:
      try:
        ll.append(float(val))
      except Exception:
        raise ConversionError, "Unable to convert float '%s'" % val
    return ll



class IntConverter(Converter):
  """Converts string to integer"""

  def convert(self, values):
    ll = []
    for val in values:
      try:
        ll.append(int(val))
      except Exception:
        raise ConversionError, "Unable to convert integer '%s'" % val
    return ll



class ComplexConverter(Converter):
  """Converts string to complex"""

  def __call__(self, strValue):
    values = strValue.split()
    if len(values) % 2:
      raise ConversionError, "Odd number of values"
    if self.nolist and len(values) != 2:
      raise ConversionError, "Too many values"
    result = self.convert(values)
    if self.nolist:
      return result[0]
    else:
      return result


  def convert(self, values):
    ll = []
    for ii in range(0, len(values), 2):
      try:
        ll.append(complex(float(values[ii]), float(values[ii+1])))
      except Exception:
        raise ConversionError, ("Unable to convert complex '(%s,%s)'"
                                % (values[ii], values[ii+1]))
    return ll



class LogicalConverter(Converter):
  """Converts string to logical"""

  def convert(self, values):
    ll = []
    for val in values:
      if val == 'T' or val == 't':
        ll.append(1)
      elif val == 'F' or val == 'f':
        ll.append(0)
      else:
        raise ConversionError, "Unable to convert logical '%s'" % val
    return ll
############################################################################
# Tagged data related objects
############################################################################

class TaggedEntry(object):
  """Represents a tagged entry with data."""

  # Converter from string for different types
  __strToValue = { "integer" : IntConverter(),
                   "real"    : FloatConverter(),
                   "complex" : ComplexConverter(),
                   "logical" : LogicalConverter()
                   }

  # Valid types
  __validTypes = __strToValue.keys()


  def __init__(self, name, type, rank, shape, strValue):
    """Instantiates an TaggedEntry object.
    name     -- name of the tagged entry
    type     -- type of the data
    rank     -- rank of the data
    shape    -- shape of the data (as tuple)
    strValue -- 
    """

    if not type in self.__validTypes:
      raise InvalidEntry(msg="Invalid data type '%s'" % type)
    self.__name = name
    self.__type = type
    self.__rank = rank
    self.__shape = shape
    try:
      self.__value = self.__strToValue[type](strValue)
    except ConversionError, msg:
      raise InvalidEntry(msg=msg)
    if shape and (len(self.__value) != reduce(lambda x,y: x*y, shape)):
      raise InvalidEntry(msg="Invalid nr. of values")


  def getName(self):
    return self.__name
  name = property(getName, None, None, "name of the entry")


  def getType(self):
    return self.__type
  type = property(getType, None, None, "type of the data in the entry")


  def getRank(self):
    return self.__rank
  rank = property(getRank, None, None, "rank of the data in the entry")


  def getShape(self):
    return self.__shape
  shape = property(getShape, None, None, "shape of the data in the entry")


  def getValue(self):
    return self.__value
  value = property(getValue, None, None, "value of the data in the entry")


  def isComparable(self, other):
    """Check if two entries are comparable"""
    return (other.name == self.name and other.type == self.type
            and other.rank == self.rank and other.shape == self.shape)




class TaggedCollection(object):
  """Contains a collection of tagged entries"""

  def __init__(self, entries):
    """file -- open file like object containing collection of tagged data"""
    self.__entryNames = []
    self.__entryLines = []
    self.__entries = []
    self.addEntries(entries)


  def addEntries(self, entries):

    for entry in entries:
      taggedLine = ":".join((entry.name, entry.type, str(entry.rank),
                             ",".join(map(str, entry.shape))))
      self.__entryNames.append(entry.name)
      self.__entryLines.append(taggedLine)
      self.__entries.append(entry)


  def getMatchingEntries(self, pattern):
    """Returns entries from the collection matching a given pattern
    pattern -- compiled regular expression
    """
    result = []
    for iEntry in range(len(self.__entries)):
      if pattern.match(self.__entryLines[iEntry]):
        result.append(self.__entries[iEntry])

    return result


  def getEntry(self, name):
    """Returns an entry with a given name from the collection
    name -- name of the entry
    """
    try:
      iEntry = self.__entryNames.index(name)
    except ValueError:
      result = None
    else:
      result = self.__entries[iEntry]

    return result


  def delEntry(self, name):
    """Deletes the specified entry from the collection
    name -- name of the entry
    """
    try:
      iEntry = self.__entryNames.index(name)
    except ValueError:
      pass
    else:
      del self.__entries[iEntry]
      del self.__entryNames[iEntry]
      del self.__entryLines[iEntry]

class ResultParser(object):
  """Parser the result files containing tagged data"""

  # Pattern for lines containing the describing tag for following data
  patTagLine = re.compile(r"""(?P<name>[^: ]+)\s*:
                              (?P<type>[^:]+):
                              (?P<rank>\d):
                              (?P<shape>(?:\d+(?:,\d+)*)*)
                              """, re.VERBOSE)


  def __init__(self, file):
    """file -- file like object containing tagged data"""
    self.__file = file


  def iterateEntries(self):
    """Generator for iterating over the entries of the data file."""

    name = None
    type = None
    rank = None
    shape = None
    value = []
    for line in self.__file.readlines():
      failed = False
      match = self.patTagLine.match(line)
      if match:

        # Append data from previous tag if present
        if name:
          try:
            yield TaggedEntry(name, type, rank, shape, " ".join(value))
          except InvalidEntry, ee:
            raise InvalidEntry(0, 0, msg=ee.msg)

        name = match.group("name")
        type = match.group("type")
        rank = int(match.group("rank"))
        if rank > 0:
          shape = tuple([ int(s) for s in match.group("shape").split(",") ])
        else:
          shape = ()
        value = []
##        iTaggedLine = iLine
      else:
        value.append(line)

    # process last entry
#    if name:
#      try:
      yield TaggedEntry(name, type, rank, shape, " ".join(value))
#      except InvalidEntry, ee:
#        raise InvalidEntry(iTaggedLine + 1, iLine, msg=ee.msg)

  entries = property(iterateEntries, None, None, "Iterator over parsed entries")

################################################################################
# End Code imported from Balint
################################################################################


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
			print "initializing noodle calculator"
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
					print 'noodle calculator: workdir "%s" does not exist, creating it.'
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
		# first input the  user-provided parameters, then specify all our own data with override
		if self.paraminclude[-4:-1].lower==".xml":
			print >> ninput, "<<! %s" % (self.paraminclude)
		else:
			print >> ninput, "<<+ %s" % (self.paraminclude)
		#keep the geometry file separate
		print >> ninput, """Geometry = GenFormat {\n <<< "input.gen" \n}"""
##		#this should hopefully give us forces and total energies at a single point
##		print >> ninput, "!Driver = ConjugateGradient{MovedAtoms=Range{1 -1} MaxSteps=0}"
		#override initial charge reuse according to options  and existence of charge file
		if self.rchr:
			if os.path.exists("charges.bin"):
				print >> ninput, "*Hamiltonian = *DFTB {!ReadInitialCharges = Yes}"
		else:
			print >> ninput, "*Hamiltonian = *DFTB {!ReadInitialCharges = No}"
		#specify system charge
		print >> ninput, "*Hamiltonian = *DFTB {!Charge = %f}" % charge
		#specify the Slater-Koster files and Max angular momenta
		sklist=[]
		mxalist=[]
		symlist,symdict=Geo.getatomsymlistdict()
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
		print >> ninput, "*Hamiltonian = *DFTB {!SlaterKosterFiles = {"
		print >> ninput, newline.join(sklist)+"}"
		print >> ninput, "!MaxAngularMomentum = {"
		print >> ninput, newline.join(mxalist)+"}}"
		#override output options to our own needs
		print >> ninput, """*Options = {
!AtomResolvedEnergies = No
!WriteResultsTag = Yes
!CalculateForces = Yes}"""
		#set pointcharges options, if specified
		if pchr:
			print >> ninput,"*Hamiltonian = *DFTB {*ElectricField ={ *PointCharges= {"
			print >> ninput,'<<< "pointcharges.xyzq"'
			print >> ninput,'}}}'
		ninput.close()



	def _prepare(self, steplabel, Geometry, charge):
		"""prepare NOODLE calculator run c.f. base class"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print "preparing noodle run"
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
		symlist,symdict=Geometry.getatomsymlistdict()
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
			print "noodle postrun cleanup and statistics"
		# first some statistics
		self.totalscf+=self.scfit
		self.totalruns+=1
		if self.verbosity>=constants.VBL_TALKY:
			print "%s: SCC iterations: %3d   ----   Total Energy: %12.6fH" % (steplabel,self.scfit,self.etot)
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
		@param : atomcount number of atoms in system (ignored!)"""
		if self.verbosity>=constants.VBL_DEBUG2:
			print "parsing noodle output"
		#overwrite atomcount, if we were working on a subgeometry
		if self.remapatoms!=None:
			realatomcount=atomcount
			atomcount=self.remapatoms[0]
		#first read the results.tag file into memory
		resultsfile=utils.compressedopen("results.tag","r")
		tagresults=TaggedCollection(ResultParser(resultsfile).entries)
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
				print "Warning: noodle SCC not converged, forces may be wrong!"
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
			print "noodle calculator shut down"

