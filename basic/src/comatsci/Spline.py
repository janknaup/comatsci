## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

##############################################################################
# Spline.py
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

from __future__ import print_function

import splext #@UnresolvedImport
import bisect
import numpy.oldnumeric as num
import math

class spline:
	"""class for 1D spline interpolation of Numeric arrays"""
	
	
	
	def __init__(self, x=None, y=None):
		"""Constructor, directly calculate coefficients, if x and y arrays are passed
		@param x: abscissa values
		@param y: ordinate values"""
		if x==None and y==None:
			self.xgrid=[]
			self.ygrid=[]
			self.x2=[]
		else:
			self.setdata(x,y)

	
	
	def setdata(self,x,y):
		"""set spline input grid and calculate input second derivatives
		@param x: abscissa values
		@param y: ordinate values"""
##		if not ((type(x) is 'array') and (type(y) is 'array') and len(x)==len(y)):
		if not len(x)==len(y):
			raise("Incompatible X and Y initialization data passed")
		elif num.sometrue(x-num.sort(x)):
			raise("x-array must be sorted")
		else:
			self.xgrid=x
			self.ygrid=y
			self.x2grid=splext.spline(x,y)



	def splint(self,x):
		"""interpolate value at x
		@param x: abscissa value"""
		if self.xgrid==[]:
			raise("No data defined")
		elif x<self.xgrid[0] or x>self.xgrid[len(self.xgrid)-1]:
			raise("interpolation value out of fitrange")
		else:
			return splext.splint(self.xgrid,self.ygrid,self.x2grid,x)
	
	
	
	def splder(self,x):
		"""return spline first derivative at x
		@param x: abscissa value"""
		if self.xgrid==[]:
			raise("No data defined")
		elif x<self.xgrid[0] or x>self.xgrid[len(self.xgrid)-1]:
			raise("interpolation value out of fitrange")
		else:
			return splext.splder(self.xgrid,self.ygrid,self.x2grid,x)
		



class vectorSpline:
	"""Class for spline interpolation of vectors
	implementation based on splext implementation but don't expect high speeds"""
	
	def __init__(self, vectorFunction=None):
		"""initialize vectoSpline
		@param vectorFunction: tuple of independent variables sequence and function value vector sequence, dimmensions of vectorFunction[0] and vectorFunction[1] must be identical, all vectorFunction[1][n] sequences must be of same len(). Independents of nodes must be sorted ascending! Default: do not store initial values (default None)		
"""
		if vectorFunction!=None:
			self.setNodes(vectorFunction)
		else:
			self._nodecount=0
			self._dimension=0



	def getDimension(self):
		return self._dimension
	dimension=property(fget=getDimension,doc="Dimension of the interpolated vectors")
	
	
	
	def getNodeCount(self):
		return self._nodecount
	nodeCount=property(fget=getNodeCount,doc="Number of nodes")
	
	
	
	def getHasNodes(self):
		"""check if necessary nodes for interpolation are present
		@return: True if >=2 nodes are stored, False otherwise"""
		if len(self._nodesx)>=2:
			return True
		else:
			return False
	hasNodes=property(fget=getHasNodes, doc="""True if nodes for interpoation are present, False otherwise""")
	
	
	
	def setNodes(self, nodes):
		"""set spline nodes for interpolation
		@param nodes: tuple of independent variables sequence and function value vector sequence, dimmensions of nodes[0] and nodes[1] must be identical, all nodes[1][n] sequences must be of same len(). Independents of nodes must be sorted ascending!"""
		# first check consitency of nodes
		if len(nodes[0])!=len(nodes[1]):
			raise ValueError,"numbers of independent variables and vectors differ"
		else:
			self._nodecount=len(nodes[0])
			self._dimension=len(nodes[1][0])
			for i in nodes[1]:
				if len(i)!=self._dimension:
					raise ValueError,"node vectors not of equal dimension"
		# check monotonicity of indepenent values
		for i in range(1,len(nodes[0])):
			if nodes[0][i]<=nodes[0][i-1]:
				raise ValueError,"node indepenent values not strictly monotonous"
		# store nodes in internal format (Doubles, one 1D-1D spline per function falue coordinate)
		self._nodesx=num.array(nodes[0],num.Float)
		temparray=num.transpose(num.array(nodes[1],num.Float))
		self._coordinateSplines=[]
		for i in range(self._dimension):
			self._coordinateSplines.append(spline(x=self._nodesx,y=temparray[i]))
		self.__arclength=None
		
		
		
	
	def splint(self, x):
		"""interpolate vector function
		@param x: independent variable value to interpolate at
		@return: interpolated vector"""
		# sanity checks
		if not self.hasNodes:
			raise "no nodes available for interpolation"
		elif x < self._nodesx[0] or x > self._nodesx[-1]:
			raise ValueError,"independent value out of node range"
		# interpolate individual coordinates
		retArray=num.zeros((self._dimension),num.Float)
		for i in range(self._dimension):
			retArray[i]=self._coordinateSplines[i].splint(x)
		return retArray
	
	
	
	def splder(self, x):
		"""interpolate vector function gradient
		@param x: independent variable value
		@return: interpolated gradient"""
		# sanity checks
		if not self.hasNodes:
			raise "no nodes available for interpolation"
		elif x < self._nodesx[0] or x > self._nodesx[-1]:
			raise ValueError,"independent value out of node range"
		# interpolate individual coordinates
		retArray=num.zeros((self._dimension),num.Float)
		for i in range(self._dimension):
			retArray[i]=self._coordinateSplines[i].splder(x)
		return retArray


	
	def getArcLength(self, tolerance=1e-6, points=None):
		"""return numerically integrated arc-length of current spline
		@param tolerance: maximum relative difference btw. arc length iterations (default 1e-6)
		@param points: 2-tuple of parameter values between which to integrate. default: integrate along whole parameter range (default None)
		@return: arc length along vector spline along specified parameter interval"""
		# sanity checks
		if not self.hasNodes:
			raise "no nodes to calculate arc length with"
		if points==None:
			a=self._nodesx[0]
			b=self._nodesx[-1]
		else:
			a=float(min(points))
			b=float(max(points))
			if a<self._nodesx[0]:
				raise ValueError,"Vectorspline: lower arc-length interval out of paramter bounds"
			if b>self._nodesx[-1]:
				raise ValueError,"Vectorspline: upper arc-length interval out of paramter bounds"
		# if total arc length is requested, we may have a stored value we can re-use
		if points==None and self.__arclength!=None:
			return self.__arclength
		#start with 1000 samples, simultaneously calculate arc length for 1000 and 500 samples
		nsamples=1000
		# calculate first sample out of order to ba able to perform
		# summation in sampling loop
		oldsample=self.splint(a)
		hiarc=0.
		loarc=0.
		tolerance=1e-6
		step=(b-a)/float(nsamples-1)
		self.arcsamples=num.zeros((nsamples,self._dimension),num.Float)
		for i in range(1,nsamples/2):
			loposition=a+(2.*i*step)
			hiposition=a+(((2.*i)-1.)*step)
			losample=self.splint(loposition)
			hisample=self.splint(hiposition)
			diff=losample-oldsample
			loarc+=num.sqrt(num.dot(diff,diff))
			diff=hisample-oldsample
			hiarc+=num.sqrt(num.dot(diff,diff))
			diff=losample-hisample
			hiarc+=num.sqrt(num.dot(diff,diff))
			oldsample=losample
			self.arcsamples[i*2]=losample
			self.arcsamples[(i*2)-1]=hisample
##		print("hi: {0:f}  lo: {1:f}".format(hiarc,loarc))
		while abs(hiarc-loarc)/hiarc>tolerance:
			nsamples*=2
			loarc=hiarc
			hiarc=0.
			step=(b-a)/float(nsamples-1)
			self.arcsamples=num.zeros((nsamples,self._dimension),num.Float)
			oldsample=self.splint(a)
			for i in range(1,nsamples):
				hiposition=a+(i*step)
				hisample=self.splint(hiposition)
				diff=hisample-oldsample
				hiarc+=num.sqrt(num.dot(diff,diff))
				oldsample=hisample
				self.arcsamples[i]=hisample
##			print("hi: {0:f}  lo: {1:f}".format(hiarc,loarc))
		# if we calculated the total arc length, store it for re-use
		if points==None:
			self.__arclength=hiarc
		return hiarc
	
	arcLength=property(getArcLength,doc="total arc length of spline at default tolerance")


class RennerSpline:
	"""Renner subspline interpolation of 2D arrays"""
	
	def __init__(self, vectors=None, closed=False):
		"""Construct Renner subspline object
		@param vectors: subscriptable of 1D Numeric arrays. Must provide at least 4 points
		@param closed: boolean specifying whether the Renner subspline curve is to be open or closed (default False)		
"""
		# initialize if user has specified vectors, otherweise don't
		if vectors==None:
			self.__initalized=False
		else:
			self.setNodes(vectors,closed)




	def reset(self):
		"""reset RennerSpline object to unitialized state"""
		self._intervalLengths=None
		self._lengthBefore=None
		self._coefficients=None
		self._tangentsLeft=None
		self._tangentsRight=None
		self._chords=None
		self._normalizedChords=None
		self._totalLength=None
		self.__initialized=False
	
	
	
	
	def setNodes(self, vectors, closed=False):
		"""calculate spline parameters from node vectors
		@param vectors: subscriptable of 1D Numeric arrays
		@param closed: boolean specifying whether the Renner subspline curve is to be open or closed (default False)		
"""
		# reset self, to be sure no old data is present
		self.reset()
		# sanity-check vectors
		if len(vectors) <5:
			raise "Renner subspline needs at least 4 points"
		# store intervals count and vector dimension
		intervals=len(vectors)-1
		dimension=len(vectors[0])
		# allocate arrays in advance
		self._intervalLengths=num.zeros(intervals,num.Float)
		self._lengthBefore=num.zeros(intervals,num.Float)
		self._coefficients=num.zeros((intervals,4,dimension),num.Float)
		self._tangentsLeft=num.zeros((intervals+1,dimension),num.Float)
		self._tangentsRight=num.zeros((intervals+1,dimension),num.Float)
		# weird storage here! We need chords s_-2,s_-1,s_0...,s_n,s_n+1
		# we store s_-2 and s_-1 at the end of the array, so we can index them as -2 and -1
		self._chords=num.zeros((intervals+4,dimension),num.Float)
		self._normalizedChords=num.zeros((intervals+4,dimension),num.Float)
		#----------------------------------------------------------------------------------
		# calculate and store chords
		# start with chords 0..n-1
		for i in range(intervals):
			self._chords[i]=vectors[i+1]-vectors[i]
		# distinguish between open and closed curve
		if closed:
			self._chords[-2]=self._chords[intervals-1]
			self._chords[-1]=self._chords[intervals]
			self._chords[intervals+1]=self._chords[0]
			self._chords[intervals+2]=self._chords[1]
		else:
			self._chords[-2]=3.*self._chords[0]-2.*self._chords[1]
			self._chords[-1]=2.*self._chords[0]-1.*self._chords[1]
			self._chords[intervals+1]=2.*self._chords[intervals]-1.*self._chords[intervals-1]
			self._chords[intervals+2]=3.*self._chords[intervals]-2.*self._chords[intervals-1]
		# calculate normalized chords
		for i in range(len(self._chords)):
			# chords can be the zero-vector, in that case normalized chord should also be zero
			chordlen=math.sqrt(num.dot(self._chords[i],self._chords[i]))
			if not chordlen<1E-12:
				self._normalizedChords[i]=self._chords[i]/chordlen
			else:
				self._normalizedChords[i]=num.zeros(dimension,num.Float)
		# chords calculation finished
		#----------------------------------------------------------------------------------
		#----------------------------------------------------------------------------------
		# calculate and store left and right tangents
		for i in range(0,intervals+1):
			alpha=math.sqrt(1-num.dot(self._normalizedChords[i-2],self._normalizedChords[i-1]))
			NE=alpha+math.sqrt(1-num.dot(self._normalizedChords[i],self._normalizedChords[i+1]))
			if NE>0:
				alpha/=NE
				self._tangentsLeft[i]=self._chords[i-1]+(alpha*(self._chords[i]-self._chords[i-1]))
				self._tangentsLeft[i]/=math.sqrt(num.dot(self._tangentsLeft[i],self._tangentsLeft[i]))
				self._tangentsRight[i]=self._tangentsLeft[i]
			else:
				self._tangentsLeft[i]=self._normalizedChords[i-1]
				self._tangentsLeft[i]=self._normalizedChords[i]
		# tangents calculation finished
		#----------------------------------------------------------------------------------
		#----------------------------------------------------------------------------------
		# calculate and store parameter interval lengths
		# lengthbefore stores the parameter length before the current interval
		for i in range(intervals):
			A=16-num.dot((self._tangentsRight[i]+self._tangentsLeft[i+1]),(self._tangentsRight[i]+self._tangentsLeft[i+1]))
			B=6.*num.dot(self._chords[i],(self._tangentsRight[i]+self._tangentsLeft[i+1]))
			C=36.*num.dot(self._chords[i],self._chords[i])
			self._intervalLengths[i]=(-B+math.sqrt(B*B+A*C))/A
			# all arrays are initialized to zeros, therefore lengthBefore[0] does not have to be set
			if i > 0:
				self._lengthBefore[i]=self._lengthBefore[i-1]+self._intervalLengths[i]
		# parameter interval lengths finished
		#----------------------------------------------------------------------------------
		#----------------------------------------------------------------------------------
		# calculate and store coefficients
		for i in range(intervals):
			#a_i
			self._coefficients[i][0]=vectors[i]
			#b_i
			self._coefficients[i][1]=self._tangentsRight[i]
			#c_i
			self._coefficients[i][2]=(3./(self._intervalLengths[i]**2))*self._chords[i]
			self._coefficients[i][2]-=(1./self._intervalLengths[i])*(2.*self._tangentsRight[i]+self._tangentsLeft[i+1])
			#d_i
			self._coefficients[i][3]=(1./self._intervalLengths[i]**2)*(self._tangentsRight[i]+self._tangentsLeft[i+1])
			self._coefficients[i][3]-=(2./self._intervalLengths[i]**3)*self._chords[i]
		# coefficients finished
		#----------------------------------------------------------------------------------
		# calculate total arc length
		self._totalLength=self._lengthBefore[intervals-1]+self._intervalLengths[intervals-1]
		# mark self as initialized
		self.__initialized=True
		#done
	
	
	
	def splint(self, t):
		"""(Renner sub-)spline interpolate vector curve at parameter value t
		@param t: parameter value at which to interpolate must be in 0..totalLength"""
		# check if we are initialized
		if not self.__initialized:
			raise "Renner Spline object not initialized before interpolation"
		# sanity-check t
		if t < 0.0 or t > self._totalLength:
			raise ValueError("interpolation point out of parameter Range")
		# find which spline interval we are in
		interval=bisect.bisect_right(self._lengthBefore,t)-1
		# project t into current interval
		t-=self._lengthBefore[interval]
		# interpolate
		S=num.zeros(num.shape(self._coefficients[0][0]),num.Float)
		for i in range(0,4):
			S+=self._coefficients[interval][i]*t**i
		return S
		# done
	
	
	
	def splder(self,t):
		"""derivative of the Renner subspline at parameter value t
		@param t: parameter value at which to derivate must be in 0..totalLength"""
		# check if we are initialized
		if not self.__initialized:
			raise "Renner Spline object not initialized before interpolation"
		# sanity-check t
		if t < 0.0 or t > self._totalLength:
			raise ValueError("interpolation point out of parameter Range")
		# find which spline interval we are in
		interval=bisect.bisect_right(self._lengthBefore,t)-1
		# project t into current interval
		t-=self._lengthBefore[interval]
		# interpolate derivative
		# dS/dt = b+1/2*c*T+1/3*d*t^2
		dS=self._coefficients[interval][1]
		dS+=1./2. * self._coefficients[interval][2]*t
		dS+=1./3. * self._coefficients[interval][3]*t*t
		return dS
		# done



	def getTotalLength(self):
		"""get the total curve length of the Renner Subspline
		@return: total length of the Subspline curve"""
		return self._totalLength
	
	totalLength=property(getTotalLength,None,None,"total curve length of the subspline")
