##############################################################################
# utils.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2012 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

import gzip,bz2,os,shutil
import time,sys,re



def compressedopen(filename, mode="r", compresslevel=0, autodetect=True):
	"""return a fileobject with specified mode. Open "filename" uncompressed if existant, otherwise try to open filename.gz or filename.bz2. If compresslevel >=1 and <=9 is specified upon open in write mode, a gzip file of specified compression is opened.
	If filename exists, compressed files are ignored, if filename does not exist, and both compressed variants exist, an exception is raised and no file is opened
	@param filename: Name of the file to be opened
	@param mode: Mode in which the file is to be opened, same as for standard file objects
	@param compresslevel: compression level for files openened for write. If >0 and <=9, a compressed file named filename.gz is opened (default 0)
	@param autodetect: Boolean, specifying, if a compressed file should be automatically searched, or not (default True)
	@return: file object, gzipFileObject or bz2FileObject
	"""
	# compression level sanity check
	if compresslevel > 9 or compresslevel <0 :
		raise ValueError("Invalid compression level specified!")
	#  check if filename already specifies a compressed file. If true, try to open that file only, using the proper compressed file object
	if filename[-3:].lower()==".gz":
		# compresslevel 0 is undefined in gzip and bz2, so work around it
		if compresslevel>0:
			return gzip.open(filename,mode,compresslevel)
		else:
			return gzip.open(filename,mode)
	elif filename[-4:].lower()==".bz2":
		if compresslevel>0:
			return bz2.BZ2File(filename,mode,compresslevel=compresslevel)
		else:
			return bz2.BZ2File(filename,mode)
	else:
		# check if an existing file is to be opened
		mode=mode.lower()
		if mode[0]=="r" or mode[0]=="a":
			# if uncompressed file exists, open that
			if os.path.exists(filename):
				return open(filename,mode)
			#otherwise check for compressed files
			else:
				bz2exist=os.path.exists(filename+".bz2")
				gzexist=os.path.exists(filename+".gz")
				# bail out if both compressed variants exist
				if bz2exist and gzexist and not autodetect:
					raise ValueError("Cannot decide wether to open {0:s} or {1:s}.".format(filename+".bz2",filename+".gz"))
				# if one of the compressed variants exists, open that:
				elif bz2exist and autodetect:
					if compresslevel>0:
						return bz2.BZ2File(filename+".bz2",mode,compresslevel=compresslevel)
					else:
						return bz2.BZ2File(filename+".bz2",mode)
				elif gzexist and autodetect:
					return gzip.open(filename+".gz",mode,compresslevel)
				# if none exists and we have mode a, create an umcompressed file for compresslevel 0 and a bzip file otherwise
				elif mode[0]=="a":
					if compresslevel==0:
						return open(filename,mode)
					else:
						return gzip.open(filename,mode,compresslevel)
				# if all else fails, raise a file not found exception
				else:
					raise ValueError("No file candiate to open for read was found for file '{0:s}'".format(filename))
		# If a file is to be opened for write, just create a normal or compressed file, based on compresslevel
		else:
			if compresslevel==0:
				return open(filename,mode)
			else:
				return gzip.open(filename+".gz",mode,compresslevel)




def uncompresscopy(source, destination):
	"""copy source file to destination. If source file is gzip or bz2 comressed, uncompress it to destination.
	<p><b>Sufficient storage space for the uncompressed file must exist at the <em>source</em> location </b></p>
	@param source: source file path to copy and uncompress, if necessary
	@param destination: destination path to copy (uncompressed) source to"""
	# just check for a known compressed format filename extension and handle accordingly
	if source[-3:].lower()==".gz":
		os.system("gzip -d {0:s}".format(source))
		shutil.copy(source[:-3],destination)
	elif source[-4:].lower()==".bz2":
		os.system("gzip -d {0:s}".format(source))
		shutil.copy(source[:-4],destination)
	# simply copy otherwise
	else:
		shutil.copy(source,destination)




def compresscopy(source,destination,compresslevel=9):
	"""compress source file and copy to destination
	@param source: source file path to compress and copy
	@param destination: destination path to copy compressed source to
	@param compresslevel: compression level to use (default 9)
	"""
	# sanity check compresslevel
	if compresslevel <1 or compresslevel > 9:
		raise ValueError("Invalid compression level specified")
	# call gzip to compress
	os.system("gzip -{0:d} {1:s}".format(compresslevel,source))
	# copy compressed file
	shutil.copy(source+".gz",destination)



#################################################################################
##  Helper classes to robustly parse noodle output
## written by Balint Aradi, copy & pasted by Jan Knaup
#################################################################################
#
#class ConversionError(Exception):
#	pass
#
#class InvalidEntry(Exception):
#	pass
#
#############################################################################
## Conversion functors
#############################################################################
#
#class Converter(object):
#	"""Base class for string to value converters"""
#
#	def __init__(self, nolist=False):
#		"""nolist -- if True, not a list, but only a single value is returned after conversion.
#		"""
#		self.nolist = nolist
#
#	def __call__(self, strValue):
#		"""strValue -- string represenation of the values to convert"""
#		values = strValue.split()
#		if self.nolist and len(values) > 1:
#			raise ConversionError, "Too many values"
#		result = self.convert(values)
#		if self.nolist:
#			return result[0]
#		else:
#			return result
#
#	def convert(self, values):
#		"""Conversion function
#		values -- list of strings representation of values to convert
#		"""
#		return values
#
#
#
#class FloatConverter(Converter):
#	"""Converts string to float"""
#
#	def convert(self, values):
#		ll = []
#		for val in values:
#			try:
#				ll.append(float(val))
#			except Exception:
#				raise ConversionError, "Unable to convert float '{0:s}'".format(val)
#		return ll
#
#
#
#class IntConverter(Converter):
#	"""Converts string to integer"""
#
#	def convert(self, values):
#		ll = []
#		for val in values:
#			try:
#				ll.append(int(val))
#			except Exception:
#				raise ConversionError, "Unable to convert integer '{0:s}'".format(val)
#		return ll
#
#
#
#class ComplexConverter(Converter):
#	"""Converts string to complex"""
#
#	def __call__(self, strValue):
#		values = strValue.split()
#		if len(values) % 2:
#			raise ConversionError, "Odd number of values"
#		if self.nolist and len(values) != 2:
#			raise ConversionError, "Too many values"
#		result = self.convert(values)
#		if self.nolist:
#			return result[0]
#		else:
#			return result
#
#
#	def convert(self, values):
#		ll = []
#		for ii in range(0, len(values), 2):
#			try:
#				ll.append(complex(float(values[ii]), float(values[ii+1])))
#			except Exception:
#				raise ConversionError, ("Unable to convert complex '({0:s},{1:s})'".format(
#					values[ii], values[ii+1]))
#		return ll
#
#
#
#class LogicalConverter(Converter):
#	"""Converts string to logical"""
#
#	def convert(self, values):
#		ll = []
#		for val in values:
#			if val == 'T' or val == 't':
#				ll.append(1)
#			elif val == 'F' or val == 'f':
#				ll.append(0)
#			else:
#				raise ConversionError, "Unable to convert logical '{0:s}'".format(val)
#		return ll
#
#############################################################################
## Tagged data related objects
#############################################################################
#
#class TaggedEntry(object):
#	"""Represents a tagged entry with data."""
#	# Converter from string for different types
#	__strToValue = { "integer" : IntConverter(),
#					"real"	: FloatConverter(),
#					"complex" : ComplexConverter(),
#					"logical" : LogicalConverter()
#					}
#	# Valid types
#
#	__validTypes = __strToValue.keys()
#
#
#	def __init__(self, name, tpe, rank, shape, strValue):
#		"""Instantiates an TaggedEntry object.
#		name	 -- name of the tagged entry
#		tpe	 -- type of the data
#		rank	 -- rank of the data
#		shape	-- shape of the data (as tuple)
#		strValue -- 
#		"""
#
#		if not tpe in self.__validTypes:
#			raise InvalidEntry(msg="Invalid data type '{0:s}'".format(tpe))
#		self.__name = name
#		self.__tpe = tpe
#		self.__rank = rank
#		self.__shape = shape
#		try:
#			self.__value = self.__strToValue[tpe](strValue)
#		except ConversionError, msg:
#			raise InvalidEntry(msg=msg)
#		if shape and (len(self.__value) != reduce(lambda x,y: x*y, shape)):
#			raise InvalidEntry(msg="Invalid nr. of values")
#
#
#	def getName(self):
#		return self.__name
#	name = property(getName, None, None, "name of the entry")
#
#
#	def getType(self):
#		return self.__tpe
#	tpe = property(getType, None, None, "tpe of the data in the entry")
#
#
#	def getRank(self):
#		return self.__rank
#	rank = property(getRank, None, None, "rank of the data in the entry")
#
#
#	def getShape(self):
#		return self.__shape
#	shape = property(getShape, None, None, "shape of the data in the entry")
#
#
#	def getValue(self):
#		return self.__value
#	value = property(getValue, None, None, "value of the data in the entry")
#
#
#	def isComparable(self, other):
#		"""Check if two entries are comparable"""
#		return (other.name == self.name and other.tpe == self.tpe
#			and other.rank == self.rank and other.shape == self.shape)
#
#
#
#
#class TaggedCollection(object):
#	"""Contains a collection of tagged entries"""
#
#	def __init__(self, entries):
#		"""file -- open file like object containing collection of tagged data"""
#		self.__entryNames = []
#		self.__entryLines = []
#		self.__entries = []
#		self.addEntries(entries)
#
#
#	def addEntries(self, entries):
#		for entry in entries:
#			taggedLine = ":".join((entry.name, entry.tpe, str(entry.rank),
#								",".join(map(str, entry.shape))))
#			self.__entryNames.append(entry.name)
#			self.__entryLines.append(taggedLine)
#			self.__entries.append(entry)
#
#
#	def getMatchingEntries(self, pattern):
#		"""Returns entries from the collection matching a given pattern
#		pattern -- compiled regular expression
#		"""
#		result = []
#		for iEntry in range(len(self.__entries)):
#			if pattern.match(self.__entryLines[iEntry]):
#				result.append(self.__entries[iEntry])
#
#		return result
#
#
#	def getEntry(self, name):
#		"""Returns an entry with a given name from the collection
#		name -- name of the entry
#		"""
#		try:
#			iEntry = self.__entryNames.index(name)
#		except ValueError:
#			result = None
#		else:
#			result = self.__entries[iEntry]
#		return result
#
#
#	def delEntry(self, name):
#		"""Deletes the specified entry from the collection
#		name -- name of the entry
#		"""
#		try:
#			iEntry = self.__entryNames.index(name)
#		except ValueError:
#			pass
#		else:
#			del self.__entries[iEntry]
#			del self.__entryNames[iEntry]
#			del self.__entryLines[iEntry]
#
#class ResultParser(object):
#	"""Parser the result files containing tagged data"""
#	# Pattern for lines containing the describing tag for following data
#	patTagLine = re.compile(r"""(?P<name>[^: ]+)\s*:
#							  (?P<type>[^:]+):
#							  (?P<rank>\d):
#							  (?P<shape>(?:\d+(?:,\d+)*)*)
#							  """, re.VERBOSE)
#
#
#	def __init__(self, infile):
#		"""file -- file like object containing tagged data"""
#		self.__file = infile
#
#
#	def iterateEntries(self):
#		"""Generator for iterating over the entries of the data file."""
#		name = None
#		tpe = None
#		rank = None
#		shape = None
#		value = []
#		for line in self.__file.readlines():
#			failed = False #@UnusedVariable
#			match = self.patTagLine.match(line)
#			if match:
#				# Append data from previous tag if present
#				if name:
#					try:
#						yield TaggedEntry(name, tpe, rank, shape, " ".join(value))
#					except InvalidEntry, ee:
#						raise InvalidEntry(0, 0, msg=ee.msg)
#
#				name = match.group("name")
#				tpe = match.group("type")
#				rank = int(match.group("rank"))
#				if rank > 0:
#					shape = tuple([ int(s) for s in match.group("shape").split(",") ])
#				else:
#					shape = ()
#				value = []
#				##iTaggedLine = iLine
#			else:
#				value.append(line)
#		# process last entry
#		if name:
#			try:
#				yield TaggedEntry(name, tpe, rank, shape, " ".join(value))
#			except InvalidEntry, ee:
#				#raise InvalidEntry(iTaggedLine + 1, iLine, msg=ee.msg)
#				raise InvalidEntry(0, 0, msg=ee.msg)
#
#	entries = property(iterateEntries, None, None, "Iterator over parsed entries")
#
#################################################################################
## End Code imported from Balint
#################################################################################


################################################################################
# Text console progress meter from some web tutorial
################################################################################

class ProgressMeter(object):
	ESC = chr(27)
	def __init__(self, **kw):
		# What time do we start tracking our progress from?
		self.timestamp = kw.get('timestamp', time.time())
		# What kind of unit are we tracking?
		self.unit = str(kw.get('unit', ''))
		# Number of units to process
		self.total = int(kw.get('total', 100))
		# Number of units already processed
		self.count = int(kw.get('count', 0))
		# Refresh rate in seconds
		self.rate_refresh = float(kw.get('rate_refresh', .5))
		# Number of ticks in meter
		self.meter_ticks = int(kw.get('ticks', 60))
		self.meter_division = float(self.total) / self.meter_ticks
		self.meter_value = int(self.count / self.meter_division)
		self.last_update = None
		self.rate_history_idx = 0
		self.rate_history_len = 10
		self.rate_history = [None] * self.rate_history_len
		self.rate_current = 0.0
		self.last_refresh = 0
		self._cursor = False
		self.reset_cursor()

	def reset_cursor(self, first=False):
		if self._cursor:
			sys.stdout.write(self.ESC + '[u')
		self._cursor = True
		sys.stdout.write(self.ESC + '[s')

	def update(self, count, **kw):
		now = time.time()
		# Caclulate rate of progress
		rate = 0.0
		# Add count to Total
		self.count += count
		self.count = min(self.count, self.total)
		if self.last_update:
			delta = now - float(self.last_update)
			if delta:
				rate = count / delta
			else:
				rate = count
			self.rate_history[self.rate_history_idx] = rate
			self.rate_history_idx += 1
			self.rate_history_idx %= self.rate_history_len
			cnt = 0
			total = 0.0
			# Average rate history
			for rate in self.rate_history:
				if rate == None:
					continue
				cnt += 1
				total += rate
			rate = total / cnt
		self.rate_current = rate
		self.last_update = now
		# Device Total by meter division
		value = int(self.count / self.meter_division)
		if value > self.meter_value:
			self.meter_value = value
		if self.last_refresh:
			if (now - self.last_refresh) > self.rate_refresh or \
				(self.count >= self.total):
					self.refresh()
		else:
			self.refresh()

	def get_meter(self, **kw):
		bar = '-' * self.meter_value
		pad = ' ' * (self.meter_ticks - self.meter_value)
		perc = (float(self.count) / self.total) * 100
		return '[{0:s}>{1:s}] {2:4.1f}%  {3:.1f}/sec'.format(bar, pad, perc, self.rate_current)

	def refresh(self, **kw):
		# Clear line
		sys.stdout.write(self.ESC + '[2K')
		self.reset_cursor()
		sys.stdout.write(self.get_meter(**kw))
		# Are we finished?
		if self.count >= self.total:
			sys.stdout.write('\n')
		sys.stdout.flush()
		# Timestamp
		self.last_refresh = time.time()


def dictionaryPrettyPrint(indict, leftskip=0):
		"""transform linking results dictionary into a string for output
		@param indict: dictionary containing linking results
		@param leftskip: optionally, prepend leftskip blank columns to the whole output (for nested dictionaries) (default 0)
		@return: formatted string of link results for screen output
		"""
		#compile labels and get maximum label length
		labels=indict.keys()
		labels.sort()
		labelsWidth=0
		for i in labels:
			if len(str(i))>labelsWidth:
				labelsWidth=len(str(i))
		# now iterate through dictionary and compile output string
		#  perverse little trick to obtain variable length left column: create a template string, setting the output width of the left string part to labelsWidth. First for the leftskip, then append format string 
		lineTemplate="".ljust(leftskip)+"{{0:{0:d}s}}: {{1:s}}".format(labelsWidth)
		#  at first compile output lines as list, join at the end
		prettyLines=[]
		for i in labels:
			# check type of datum to key i. use special output format for dictionaries
			if isinstance(indict[i],dict):
				content="\n"+dictionaryPrettyPrint(indict[i],labelsWidth+3)
			else:
				content=str(indict[i])
			prettyLines.append(lineTemplate.format(i,content))
		# compile output lines to string and return
		return "\n".join(prettyLines)
