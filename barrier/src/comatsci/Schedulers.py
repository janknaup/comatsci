##############################################################################
# Schedulers.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup <Knaup@bccms.uni-bremen.de>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

SCHEDSTATUS_READY=0
SCHEDSTATUS_BUSY=1
SCHEDSTATUS_JOBERROR=2
SCHEDSTATUS_UNPREPARED=4
SCHEDSTATUS_DEAD=5

SCHEDSTATUSDICT={
	SCHEDSTATUS_READY:	"ready",
	SCHEDSTATUS_BUSY:	"busy",
	SCHEDSTATUS_JOBERROR:	"job error",
	SCHEDSTATUS_UNPREPARED:	"unprepared",
	SCHEDSTATUS_DEAD: "shut down",
}

#import verbosity levels from outside, if available, otherwise use own definitions
try:
	import constants
except:
	class constants:
		"""default Verbosity level definitions"""
		# Verbosity levels:
		# only fatal errors
		VBL_SILENCE=-1
		# fatal errors and warnings
		VBL_QUIET=0
		# VBL_QUIET+status and progress reports
		VBL_NORMAL=1
		#VBL_NORMAL+extended status reports
		VBL_TALKY=2
		#VBL_TALKY+basic debug information
		VBL_DEBUG1=3
		#VBL_DEBUG1+extended debug information
		VBL_DEBUG2=4


import os, sys
import time
import copy

class Scheduler:
	"""Base class to derive calculation scheduler subclasses from."""
	
	
	def __init__(self, verbosity=None):
		"""Initialize internal bookkeeping and statistics counters etc.
		(Base class constructor should be called first in derived class constructors.)
		@param verbosity:  Verbosity level, if none choose VBL_NORMAL  (default None)		
"""
		self._status=SCHEDSTATUS_UNPREPARED
		self._jobcounter=0
		self._iterationcounter=0
		self._walltimer=0.0
		self._cputimer=0.0
		if verbosity!=None:
			self._verbosity=verbosity
		else:
			self._verbosity=constants.VBL_NORMAL
	
	
	
	def perform(self, schedule,progressfunction=None):
		"""Perform calculations for the given list of input data, return a list of result dictionaries.
		@param schedule:  list of input data to pass to single workers 
		@param: progressfunction=None  lFor progress reporting purposes. If ! progressfunction is called after each job on schedule is finished with an argument indicating the percentage of jobs finished.  (default None)		
"""
		#check if scheduler is iterationCountready
		if self._status!=SCHEDSTATUS_READY:
			#if not, report fatal error
			raise("fatal scheduler error: scheduler %s"%SCHEDSTATUSDICT[self._status])
		else:
			if self._verbosity>=constants.VBL_TALKY:
				print "Executing a schedule of %d jobs." % (len(schedule))
			self._status=SCHEDSTATUS_BUSY
			results=self._realperform(schedule, progressfunction)
			self._status=SCHEDSTATUS_READY
			if self._verbosity>=constants.VBL_TALKY:
				print "Finished schedule."
			return results
	
	
	
	def _realperform(self, schedule, progressfunction):
		"""Really perform calculations for the given list of input data, return a list of result dictionaries. Must be reimplemented by subclass
		@param schedule:  list of input data to pass to single workers 
		@param: progressfunction=None  lFor progress reporting purposes. If ! progressfunction is called after each job on schedule is finished with an argument indicating the percentage of jobs finished.  (default None)		
"""
		pass
	
	

	def getJobCount(self):
		"Return total number of jobs finished by scheduler"
		return self._jobcounter
	jobCount=property(fget=getJobCount,doc="Total number of jobs finished by scheduler")
	
	
	
	def getIterationCount(self):
		"Retrun total number of iterations performed under control of scheduler"
		return self._iterationcounter
	iterationCount=property(getIterationCount,doc="Total number of iterations performed under control of scheduler")
	
	
	
	def getWallTime(self):
		"Retrun total wall time spent under control of scheduler"
		return self._walltimer
	wallTime=property(getWallTime,doc="Total wall time spent under control of scheduler")
	
	
	
	def getCPUTime(self):
		"Retrun total CPU time spent under control of scheduler"
		return self._cputimer
	cpuTime=property(getCPUTime, doc="Total CPU time spent under control of scheduler")
	
	
	
	def getStatus(self):
		"Retrun status of scheduler"
		return self._status
	status=property(getStatus,doc="status of scheduler")
	
	
	
	def shutdown(self):
		"""shutdown the scheduler and all associated workers"""
		#call the worker's shutdown method, ignore if not present
		try:
			self._worker.shutdown()
		except:
			pass
		self._status=SCHEDSTATUS_DEAD
		if self._verbosity>=constants.VBL_DEBUG1:
			print "Scheduler shut down"



class serialScheduler(Scheduler):
	"""Scheduler that performs jobs on the schedule one after another"""
	
	import time
	
	
	
	def __init__(self, worker, verbosity=None):
		"""Initialize serial scheduler
		@param worker:  Class that performs the work and returns the 
		results. Worker must have a run method that accepts input data sets 
		and a getresults method that returns a dictionary of results. Worker 
		should have an iterations property for iterations statistics to work.
		A shutdown method of the worker is called upon scheduler shutdown, 
		if present. 
		@param verbosity:  Verbosity level, if none choose  (default None)
		VBL_NORMAL"""
		Scheduler.__init__(self,verbosity)
		self._worker=worker
		self._status=SCHEDSTATUS_READY
	
	
	
	def _realperform(self, schedule, progressfunction):
		"""perform one job from schedule after another
		@param schedule:  list of input data to pass to single workers 
		@param: progressfunction=None  lFor progress reporting purposes. If ! progressfunction is called after each job on schedule is finished with an argument indicating the percentage of jobs finished.  (default None)		
"""
		#initialize results container
		results=[]
		#for progress reporting
		totaljobs=len(schedule)
		#now just do the work
		for i in range(len(schedule)):
			if self._verbosity>=constants.VBL_DEBUG1:
				print "Serial scheduler: starting job %d of %d" %(i, totaljobs)
			#measure CPU-time inside walltime bracket
			startwalltime=time.time()
			startcputime=time.clock()
			self._worker.run(schedule[i])
			cputime=time.clock()-startcputime
			walltime=time.time()-startwalltime
			self._cputimer+=cputime
			self._walltimer+=walltime
			#progress reporting
			if progressfunction!=None:
				progressfunction(int(i*100/totaljobs))
			if self._verbosity>=constants.VBL_DEBUG1:
				print "Serial scheduler: finished job after %10.4f seconds CPU time (%10.4f seconds wall time)." %(
					cputime, walltime)
			#worker may not implement an iteration count property
			try:
				self._iterationcounter+=self._worker.iterations
			except:
				if self._verbosity>=constants.VBL_DEBUG2:
					print "Serial scheduler: could not get iteration count from worker, assuming 1 iteration"
				self._iterationcounter+=1
			self._jobcounter+=1
			results.append(self._worker.getresults())
		return results





class mpiScheduler(Scheduler):
	"""Distribute jobs across several nodes via MPI. Relies on pypar package. All MPI initialization and communication is done transparently."""
	
	INITTAG   =   0
	WORKTAG   =   1
	RESULTTAG =   2
	DIETAG    = 254

	
	def __init__(self, initializer, initargument, Verbosity=None):
		"""Initialize MPI scheduler
		@param initializer:  A function that returns an instance to perform the work and returns the results. Worker must have a run method that accepts input data sets and a getresults method that returns a dictionary of results. Worker should have an iterations property for iterations statistics to work. 
		@param initoption:  argument that is passed to initializer when initializing the worker. Must be deep-copyable. 
		@param Verbosity:  Verbosity level, if none choose VBL_NORMAL  (default None)		
"""
		#Base class constructor
		Scheduler.__init__(self,Verbosity)
		#store the initializer in case we are slave
		self.__initializer=initializer
		#now initialize MPI
		if self._verbosity >= constants.VBL_DEBUG1:
			print "MPI scheduler: initializing MPI."
		try:
			self.pypar = __import__("pypar")
		except:
			raise 'Module pypar must be present to run parallel'
		self.__MPI_numproc=self.pypar.size()
		self.__MPI_myid =   self.pypar.rank()
		self.__MPI_node =    self.pypar.Get_processor_name()
		if self.__MPI_myid==0 and self.__MPI_numproc==1:
			print "MPI scheduler: No slaves present, no work can be performed"
			raise "no MPI slaves"
		elif self.__MPI_myid==0 and self.__MPI_numproc==2 and self._verbosity >= constants.VBL_NORMAL:
			print "*\n*\n*\nMPI scheduler: *** Warning, only one slave is present, no parallel execution possible! ***\n*\n*\n*"
		if self.__MPI_myid==0 and self._verbosity >= constants.VBL_DEBUG1:
			print "MPI scheduler: MPI initialized successfully"
		if self._verbosity >= constants.VBL_TALKY:
			print "MPI scheduler: %s: I am node %d of %d" % (self.__MPI_node, self.__MPI_myid, self.__MPI_numproc)
		if self._verbosity >=constants.VBL_NORMAL:
			if self.__MPI_myid==0:
				print "MPI scheduler: %s: I am master" % (self.__MPI_node)
			elif self._verbosity >=constants.VBL_TALKY:
				print "MPI scheduler: %s: I am slave" % (self.__MPI_node)
		# If we are master, initialize the slaves and return, otherwise switch to slave mode
		if self.__MPI_myid==0:
			if self._verbosity >=constants.VBL_NORMAL:
				print "MPI scheduler: initializing slaves"
			initmsg={
				"rundir"		: os.path.abspath('.'),
				"initarg"    : copy.deepcopy(initargument),
				}
			for i in range(1, self.__MPI_numproc):
				self.pypar.send(initmsg, i, tag=self.INITTAG)
			if self._verbosity >=constants.VBL_NORMAL:
				print "MPI scheduler: Sent init messages to slaves."
			#only the master controls the scheduler's status
			self._status=SCHEDSTATUS_READY
		else:
			self.__slave()




	def __slave(self):
		"""MPI slave"""
		# wait for init message from master and initialize
		initmsg, status = self.pypar.receive(0, tag=self.INITTAG, return_status=True)
		if self._verbosity >=constants.VBL_DEBUG2:
			sys.stderr.write("[SLAVE %d]: received init message: %s"%(self.__MPI_myid,initmsg))
		os.chdir(initmsg["rundir"])
		#create the local worker instance
		self._worker=self.__initializer(initmsg["initarg"])
		#now, be a good slave
		#i.e.: infinitelty wait for work assignments and work them off
		while 1:
			work=None
			work, status = self.pypar.receive(0, tag=self.pypar.any_tag, return_status=True) 
			if self._verbosity>=constants.VBL_DEBUG2:
				sys.stderr.write("[SLAVE %d]: received work '%s' with tag '%d' from node '%d'\n"\
					%(self.__MPI_myid, work, status.tag, status.source))
			#gracefully shutdown this slave, if DIE message is received
			if (status.tag == self.DIETAG):
				if self._verbosity >= constants.VBL_DEBUG1:
					sys.stderr.write("[SLAVE %d]: received termination from node '%d'\n" %(MPI_myid, 0))
				#if worker has a shutdown method, call it, otherwise ignore
				try:
					self._worker.shutdown()
				except:
					pass
				self.pypar.Finalize()
				if self._verbosity>=constants.VBL_TALKY:
					print "[SLAVE %d]: MPI environment finalized."%(self.__MPI_myid)
				sys.exit(0)
			#now for the work part
			elif (status.tag == self.WORKTAG):
				startwalltime=time.time()
				startcputime=time.clock()
				self._worker.run(work["work"])
				cputime=time.clock()-startcputime
				walltime=time.time()-startwalltime
				self._cputimer+=cputime
				self._walltimer+=walltime
				if self._verbosity>=constants.VBL_DEBUG1:
					print "[SLAVE %d]: finished job after %10.4f seconds CPU time (%10.4f second wall time)." %(
						self.__MPI_myid,cputime, walltime)
				result={"workresult":self._worker.getresults(),
					"cputime": cputime,
					"walltime": walltime}
				#see if the worker has an iterations property, otherwise default to 1 iteration
				try:
					result["iterations"]=self._worker.iterations
				except:
					result["iterations"]=1
				self.pypar.send(result, 0, tag=self.RESULTTAG)
				if self._verbosity >=constants.VBL_DEBUG2:
					sys.stderr.write("[SLAVE %d]: sent result '%s' to node '%d'\n" %(self.__MPI_myid, result, 0))
			#an unexpected tag is, of course, a fatal error
			else:
				sys.stderr.write("[SLAVE %d]: unexpected tag received, aborting" %(self.__MPI_myid))
				self.pypar.Finalize()
				sys.exit(1)



	def getNumProc():
		return self.__MPI_numproc
	numProc=property(getNumProc,doc="Number of MPI processors available to scheduler")
	
	
	
	def _realperform(self, schedule, progressfunction):
		"""distribute jobs to the slaves. Slaves and their workers
		must have been initialized before.
		@param schedule:  list of input data to pass to single workers 
		@param: progressfunction=None  lFor progress reporting purposes. If ! progressfunction is called after each job on schedule is finished with an argument indicating the percentage of jobs finished.  (default None)		
"""
		# make a copy of the schedule, as we work destructively on itemgetter
		myschedule=copy.deepcopy(schedule)
		# dictionary to keep track which job went to which slaves
		# necessary because we cannot expect the slaves to finish in order
		whatswhere={}
		# we will store the results in a dictionary first and then sort them into
		# a list before we return them
		resultsdict={}
		# bookkeeping
		numslaves=self.__MPI_numproc-1
		numjobs=len(myschedule)
		finished=0
		# first send out jobs to all slaves, if enough work is available
		for i in range(min(numjobs,numslaves)):
			#store which job is being sent to which slave
			#this is necessary, because jobs need not finish in order
			what=len(myschedule)-1
			whatswhere[i+1]=what
			#now get the job fromn the list and send it
			job=myschedule.pop()
			work = { "work": job}
			self.pypar.send(work, i+1, tag=self.WORKTAG) 
			if self._verbosity >= constants.VBL_DEBUG1:
				sys.stderr.write("[MASTER]: sent work '%s' to node '%d'\n" %(work, i+1))
		# inital batch of work sent out, now process remaining jobs
		while(len(myschedule)>0):
			# wait for results from slaves
			status = self.__receiveresult(resultsdict, whatswhere)
			#report progress
			if progressfunction!=None:
				progressfunction(int(finished*100/numjobs))
			what=len(myschedule)-1
			whatswhere[status.source]=len(myschedule)-1
			job=myschedule.pop()
			work = { "work": job}
			self.pypar.send(work, status.source, tag=self.WORKTAG)
			if self._verbosity >= constants.VBL_DEBUG1:
				sys.stderr.write("[MASTER]: sent job '%s' to node '%d'\n" %(work, status.source))
		# now that all work has been sent out, we must collect the remaining results
		while(len(whatswhere)>0): 
			self.__receiveresult(resultsdict, whatswhere)
			#report progress
			if progressfunction!=None:
				progressfunction(int(finished*100/numjobs))
		# finally, we must sort our results and return them as a neat list
		retresults=[]
		for i in range(numjobs):
			retresults.append(resultsdict[i])
		return retresults



	def __receiveresult(self, resultsdict, whatswhere):
		"""receive a result from the slaves, store the statistics information, store the result in resutlsdict and clean up whatswhere. return MPI status dict
		@param resultsdict:  The jobindex:result dictionary in which to store the received result 
		@param whatswhere:  The slaveid:jobindex dictionary from which to extract the index of the returned job """
		result=None
		result, status = self.pypar.receive( self.pypar.any_source, tag=self.RESULTTAG, return_status=True) 
		# a message at this point means, a calculation has finished
		self._jobcounter+=1
		self._iterationcounter+=result["iterations"]
		self._cputimer+=result["cputime"]
		self._walltimer+=result["walltime"]
		if self._verbosity > constants.VBL_DEBUG2:
			sys.stderr.write("[MASTER]: received result '%s' from node '%d'\n" %(result["workresult"], status.source))
		resultsdict[whatswhere[status.source]]=result["workresult"]
		del whatswhere[status.source]
		return status



	def shutdown(self):
		"""send shutdown messages to all slaves and shutdown local worker"""
		#send shutdown messages to slaves
		if self._verbosity >=constants.VBL_NORMAL:
			print "MPI scheduler: shutting down slaves"
		stopmsg={
			"rundir"		: os.path.abspath('.'),
			}
		for i in range(1, self.__MPI_numproc):
			self.pypar.send(stopmsg, i, tag=self.DIETAG)
		if self._verbosity >=constants.VBL_NORMAL:
			print "MPI scheduler: Sent shutdown messages to slaves."
		#shutdown local worker, ignore if no shutdown method
		try:
			self._worker.shutdown()
		except:
			pass
		#call the base class shutdown
		self.pypar.Finalize()
		if self._verbosity>=constants.VBL_NORMAL:
			print "MPI scheduler: [master] MPI environment finalized"
		Scheduler.shutdown(self)






class threadScheduler(Scheduler):
	"""Do the work parallel on the local using posix threads"""

	
	def __init__(self, initializer, initargument, maxworkthreads=0, Verbosity=None):
		"""Initialize MPI scheduler
		@param initializer:  A function that returns an instance to perform the work and returns the results. Worker must have a run method that accepts input data sets and a getresults method that returns a dictionary of results. Worker should have an iterations property for iterations statistics to work. 
		@param initoption:  argument that is passed to initializer when initializing the worker. Must be deep-copyable. 
		@param maxworkthreads:  Maximum number of concurrent worker threads to start. No limit of 0.  (default 0)
		@param Verbosity:  Verbosity level, if none choose VBL_NORMAL  (default None)		
"""
		#Base class constructor
		Scheduler.__init__(self,Verbosity)
		#store the initializer and initargument
		self.__initializer=initializer
		self.__initargument=initargument
		#store maximum thread count
		self.__maxworkthreads=maxworkthreads
		self.threading = __import__("threading")
		# we need some mutexes here:
		#  reentrant lock for access on scheduler statistics data
		self.__statslock=self.threading.RLock()
		#  reentrant lock for access on results data
		self.__resultslock=self.threading.RLock()
		#  event to signal finishing of one thread
		self.__finishedevent=self.threading.Event()
		#  semaphore to ensure that only one thread at a time signals finished
		#  (only necessary in case of very short worke runtime, but nonetheless...
		self.__finishedlock=self.threading.Semaphore(1)
		if self._verbosity>=constants.VBL_DEBUG1:
			print "threading scheduler: mutlithreading initialized"
		# now we're ready
		self._status=SCHEDSTATUS_READY
		return



	def _realperform(self, schedule, progressfunction=None):
		"""really perform the schedule in multithreading mode
		arguments: c.f. base class"""
		#first determine the number of threads to initially invoke
		if self.__maxworkthreads==0 or len(schedule)<self.__maxworkthreads:
			maxworkers=len(schedule)
		else:
			maxworkers=self.__maxworkthreads
		#initialize bookkeeping
		jobswaiting=len(schedule)
		jobsfinished=0
		numjobs=jobswaiting #avoid one call to len :-P
		#make a dictionary for storing the results
		resultsdict={}
		#initialize storage for the thread objects
		workthreads=[]
		#do not allow worker threads to signalize finished,
		#before we are ready for them
		self.__finishedlock.acquire()
		#now start the initial threads
		# we measure wall time outside the threads
		startwalltime=time.time()
		for i in range(maxworkers):
			if self._verbosity>=constants.VBL_DEBUG1:
				print "threading scheduler: starting job %d of %d" %(i+1, numjobs)
			#create thread object and immideately start the thread
			workthreads.append(self.threading.Thread(target=self.__threadworker,
				args=(i,resultsdict,schedule[i])))
			workthreads[-1].start()
			jobswaiting-=1
		#now wait for threads to finish and replace finished threads, 
		#as long as there is still work left to do_long
		self.__finishedlock.release() #allow the first thread to signalized finished
		while jobsfinished<numjobs:
			self.__finishedevent.wait()
			jobsfinished+=1
			# if there are still jobs pending, start a new thread replacing the
			# just finished one
			if jobswaiting>0:
				currentjob=numjobs-jobswaiting
				if self._verbosity>=constants.VBL_DEBUG1:
					print "threading scheduler: starting job %d of %d" %(currentjob+1, numjobs)
				workthreads.append(self.threading.Thread(target=self.__threadworker,
					args=(currentjob,resultsdict,schedule[currentjob])))
				workthreads[-1].start()
				jobswaiting-=1
			# now clear the finished signal and prepare to accept the next one
			self.__finishedevent.clear()
			self.__finishedlock.release()
		#store elapsed walltime
		self._walltimer+=time.time()-startwalltime
		#sort the results and return
		retresults=[]
		for i in range(len(schedule)):
			retresults.append(resultsdict[i])
		return retresults
	
	
	
	def getMaxWorkThreads(self):
		"""return the maximum number of worker threads"""
		return self.__maxworkthreads
	maxWorkThreads=property(getMaxWorkThreads,doc="Maximum number of concurrent worker threads to invoke")



	def __threadworker(self, jobid, resultsdict, jobarguments):
		"""helper function executed inside each single thread"""
		if self._verbosity>=constants.VBL_DEBUG2:
			print "threading scheduler: starting worker thread with input data %s" %(jobarguments)
		#initialize thread local data
		mydata = self.threading.local()
		mydata.startcputime=time.clock()
		#we need our own instance of the worker here
		mydata.myworker=self.__initializer(self.__initargument)
		#now do the actual work
		mydata.results=mydata.myworker.run(jobarguments)
		mydata.cputime=time.clock()-mydata.startcputime
		#store results with mutex
		self.__resultslock.acquire()
		resultsdict[jobid]=mydata.myworker.getresults()
		self.__resultslock.release()
		#store statistics info with mutex
		self.__statslock.acquire()
		# store iterations if provided by worker, otherwise add 1 iteration
		try:
			self._iterationcounter+=mydata.myworker.iterations
		except:
			if self._verbosity>=constants.VBL_DEBUG2:
				print "threading scheduler: could not get iteration count from worker, assuming 1 iteration"
			self._iterationcounter+=1
		# increase various statistics counters
		self._cputimer+=mydata.cputime
		self._jobcounter+=1
		self.__statslock.release()
		# signalize that we have finished, make sure that only one thread at a
		# time reports finished by acquiring semaphore
		self.__finishedlock.acquire()
		self.__finishedevent.set()
		return




class testworker:
	"""dummy worker class to test schedulers, sleeps for 1 second and returns the input as result"""
	
	def __init__(self):
		pass
	
	
	def run(self, work):
		time.sleep(1)
		self.result=work
	
	
	def getresults(self):
		return self.result
	
	
	def getiterations(self):
		return 2
	iterations=property(getiterations)
