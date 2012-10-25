## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

#
# Growstring.py
#
# Class inplementing a growing string approximation to a reaction path
# Author: Jan M. Knaup <knaup@bccms.uni.bremen.de>
#

import numpy.oldnumeric as num
from . import path,constants,spline

import bisect


class linearInterpolator():
	"""class providing linear interpolation within a set of vectors"""
	
	def __init__(self,vectors=None):
		"""initialize linear vector interpolator
		@param vectors: initial set of vectors (default None)		
"""
		if vectors!=None:
			self.setNodes(vectors)
		else:
			self.__isInitialized=False
	
	
	
	def setNodes(self,vectors):
		"""store Vector nodes and prepare interpolation"""
		# initizalize internal data storage
		self._dimension=num.shape(vectors[0])
		# linear equation is y_i=a_i*x+b_i
		self._b=num.array(vectors[:-1],num.Float)
		self._a=num.zeros((len(vectors)-1,self._dimension),num.Float)
		# we need to store the total length of the cure before each segment, to find the right segment for interpolation later
		# lengthbefore[-1] stores the total curve length
		self._lengthBefore=num.zeros((len(vectors),),num.Float)
		# iterate through vectors and store curve parameters
		for i in range(len(vectors)-1):
			segment=vectors[i+1]-vectors[i]
			length=num.sqrt(num.dot(segment,segment))
			self._a[i+1]=segment/length
			self._lengthBefore[i+1]=self._lengthBefore[i]+length
		# finally, declare ourselves as initialized
		self.__isInitialized=True
	
	
	
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------	
	def getTotalLength(self):
		"""get the total curve length
		@return: total curve length
		"""
		# total length is stored in _lengthBefore Array, cf. self.detNodes()
		if self.__isInitialized:
			return self._lengthBefore[-1]
		else:
			raise(AttributeError,"Trying to get total length of an uninitialized linear curve")
	
	totalLength=property(getTotalLength,None,None,"total curve length of the piecewise linear curve")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------



	def splint(self,x):
		"""return interpolated vector at length x - function name chosen to be api-compatible with Spline interpolators
		@param x: Flating point variable [0..curve length]
		@return: interpolated vector at lenght x along piecewise linear curve
		"""
		# generate error if called on unitialized interpolator
		if not self.__isInitialized:
			raise(AttributeError,"attempt to interpolate an uninitialized linear curve.")
		# sanity-check x to fall within parametrized range
		if x < 0:
			raise(ValueError,"negative curve-length value or parametrization specified.")
		elif x > self.totalLength:
			raise(ValueError,"interpolation point out of parameter range.")
		# find segment we are in by bisection
		pos=bisect.bisect_left(self._lengthBefore,x)
		# now calculate and returnt the interpolated vector from the appropriate inteval parameters
		return self._a[pos]*(x-self._lengthBefore[pos])+self._b[pos]



	def splder(self,x):
		"""look up first derivative of piecewise linear curve - named for api compatibility with spline interpolators
		@param x: curve length position at which to derivate
		@return: derivativ of piecewise linear curve at curve length position x
		"""
		# generate error if called on unitialized interpolator
		if not self.__isInitialized:
			raise(AttributeError,"attempt to interpolate an uninitialized linear curve.")
		# sanity-check x to fall within parametrized range
		if x < 0:
			raise(ValueError,"negative curve-length value or parametrization specified.")
		elif x > self.totalLength:
			raise(ValueError,"interpolation point out of parameter range.")
		# find segment we are in by bisection
		pos=bisect.bisect_left(self._lengthBefore,x)
		# the derivative is stored as parameter a, so just return the proper intevral value
		return self._a[pos]
	
	
	
	def reset(self):
		"""reset the linear interpolator to uninitialized state for reuse"""
		del self._dimension
		self._b=None
		self._a=None
		self._lengthBefore=None
		self.__isInitialized=False




class growingString(path.Reactionpath):
	"""A growing string representation of a chemical transition
	"""
	
	def __init__(self, checkpointdir='checkpoint', fixedatoms=None,cmode='d',fmax=0.01,ifrms=0.01,imaxit=0,
				charge=0.0,verbosity=constants.VBL_SILENCE,interpolationMode='a',targetDistance=1.0,growForce=0.1,symmetricGrowth=True,continuityTolerance=1.1):
		"""construct a growing string representation object
		@param: checkpointdir='checkpoint' directory name for checkpoint storage
		@param fixedatoms: list of atom indices to keep fixed (default None)
		@param cmode: energies and forces calculation mode
			<ul>
			<li>s serial</li>
			<li>p parallel</li>
			<li>d external scheduler</li>
			</ul>
		@param fmax: max normal force convergence criterion (default 0.01)
		@param frms: rms normal force convergence criterion (default 0.01)
		@param maxit: maximum number of iterations to perform (default 0)
		@param charge: system charge in electrons (default 0.0)
		@param: verbosity=VBL_SILENCE Verbosity level. Absolute silence by default
		@param: interpolationMode='a' interpolation mode for path parametrization and growth, can be
			<ul>
			<li>a automatic - start with linear interpolation, switch to Renner Subspline as soon as minimum number of Nodes is reached</li>
			<li>l linear - always linarly interpolate between nodes </li>
			<li>r Renner Subspline - use Renner Subspline to interpolate. Requires at least 5 geometries</li>
			</ul>
		@param targetDistance: (maximum) length of string between nodes until which to grow string (default 1.0)
		@param growForce: maximum normal fore on the end geometries of discontinuous string, before growth (default 0.1)
		@param symmetricGrowth: Ff true, grow new intermediate configuration at both discontinuous ends after both meed growForce. If False, grow new image immediately, when open any end meets growForce (default True)
		@param continuityTolerance: maximum length of a single string segment length to accept for continuous string with low parametrization density. Generate error if any segment length is > (average length)*continuityTolerance. (default 1.1)
		"""
		# Reactionpath demands an interable of fixed atom indices
		if fixedatoms==None:
			fixedatoms=[]
		# call base class constructor
		Path.Reactionpath.__int__(self,checkpointdir,fixedatoms,cmode,fmax,ifrms,imaxit,charge,verbosity) #@UndefinedVariable
		# store parameters to attributes
		self.interpolationMode=interpolationMode
		self.targetDistance=targetDistance
		self.growForce=growForce
		self.symmetricGrowth=symmetricGrowth
		self.continuityTolerance=continuityTolerance
		# initialize geometries interpolator
		# we construct an empty string here, so construct a linear interpolator in mode 'a'
		if interpolationMode in ('a','l'):
			self.__interpolator=linearInterpolator()
		elif interpolationMode=='r':
			self.__interpolator=Spline.RennerSpline()
		else:
			raise(ValueError,"Invalid growing string interpolation mode specified.")



	def checkContinuity(self):
		"""check if string is continuous. Continuity is defined as max(segmentLengths) < (average length)*continuityTolerance
		@return: True if string is continuous"""
