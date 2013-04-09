#!/usr/bin/python

import time,signal,random

import os,sys

versionname="python%d.%d" % sys.version_info[0:2]

sys.path.append(os.path.join(os.path.dirname(sys.argv[0]),
		"../lib/",versionname,"site-packages"))

import pypath

from optparse import OptionParser

import ConfigParser

# Catch SIGINT and SIGXCPU to gracefully shut down

stopsignal=None
startdir=os.path.realpath(".")

def stophandler(signum, frame):
	"""Catch Signals at which we want to stop gracefully, but abort semi hard, if signal is encountered again"""
	global stopsignal
	global startpath
	if stopsignal==None:
		if signum==signal.SIGINT:
			stopsignal="SIGINT"
		elif signum==signal.SIGXCPU:
			stopsignal="SIGXCPU"
		else:
			stopsignal=="UNKNOWN"
		print "Signal %s encountered, will stop after this NEB iteration" %(stopsignal)
		curdir=os.path.realpath(".")
		os.chdir(startpath)
		sf=open("STOP_PYPATH","w")
		sf.write(stopsignal)
		sf.close()
		os.chdir(curdir)
	else:
		print "Okay, Okay, will shutdown right now."
		sys.exit(signum)



def calcinitializer(options):
	"""initialize calculator"""
	#initialize calculator
	if options.calculator=="dftb":
		calc=pypath.Calculators.dftbcalc(optionfname="pypath.ini",verbosity=options.verbosity)
	elif options.calculator=="siesta":
		calc=pypath.Calculators.siestacalc(optionfname="pypath.ini",verbosity=options.verbosity)
	elif options.calculator=="noodle":
		calc=pypath.Calculators.noodlecalc(optionfname="pypath.ini",verbosity=options.verbosity)
	elif options.calculator=="gaussian":
		calc=pypath.Calculators.gaussiancalc(optionfname="pypath.ini",verbosity=options.verbosity)
	elif options.calculator=="muellerbrown":
		calc=pypath.Calculators.muellerBrownCalc(verbosity=options.verbosity)
	else:
		raise Exception("unknown calculator")
	return calc



def nebrun(mode="s"):
	"""run the NEB path search"""

	options, config, args = readconfig_options(mode="s")

	# Be Yakky
	
	SCHEDNAMES={
		's'	:	"serial",
		'p' : "MPI",
	}
	
	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print """
                 ******************************
                 *** pastafarian %10s ***
                 ******************************

PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis

using the following options:
     calculator           : %s
     scheduler            : %s
     maximim interations  : %d
     force tol.           : %e
     rms force tol.       : %e
     step width           : %e
     system charge        : %e\n\n""" % (pypath.constants.VERSION,options.calculator,
			SCHEDNAMES[options.scheduler],options.maxit,options.forcetol,
		options.rmstol,options.stepwidth,options.charge)

	if options.fixedatoms!="":
		tmp=options.fixedatoms.split(',')
		fixatms=[ int(s) for s in tmp ]
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "The following atoms will be kept fixed:"
			print fixatms
	else:
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print " *** WARNING: No atoms will be kept fixed. This is usually a bad idea"
		fixatms=[]

	#inititalize path and read geometry files
	fmode=config.get("NEB","fmode")
	tangmode=config.get("NEB","tangmode")
	climber=config.getint("NEB","climber")
	springk=config.getfloat("NEB","springk")
	relmode=config.get("NEB","relmode")



	
	#initialize the scheduler
	if options.scheduler.lower()=="s":
		scheduler=pypath.Schedulers.serialScheduler(calcinitializer(options), options.verbosity)
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "initialized serial scheduler"
	elif options.scheduler.lower()=="p":
		scheduler=pypath.Schedulers.mpiScheduler(calcinitializer,options, options.verbosity)
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "initialized parallel scheduler."
	elif options.scheduler.lower()=="t":
		scheduler=pypath.Schedulers.threadScheduler(calcinitializer, options, options.maxthreads, options.verbosity)
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "inialized threading schduler with maximum %d concurrent threads" % scheduler.maxWorkThreads
	else:
		raise "unknown scheduler"


	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print "%s calculator initialized" % options.calculator
		

	if options.verbosity >= pypath.constants.VBL_TALKY:
		print "initializing NEB path"
	path=pypath.NEBPath.NEBPath("checkpoint",fixatms,options.stepwidth,fmode,tangmode,
	climber,springk,'d',options.forcetol,options.rmstol,options.maxit,relmode,options.charge,options.verbosity)

	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print "reading input geometries"
	#readconfig_options() has already checked that we have either >=3 single geometry files
	#or one .fmg file, so one argument means one multi-geometry .fmg file
	if len(args)!=1:
		for i in args:
			path.appendgeofile(i)
			if options.verbosity >= pypath.constants.VBL_NORMAL:
				print i
	else:
		path.readfmgpath(args[0])
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "read fmg path from %s" % args[0]
	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print "read %d geometries" % (path.numimages())


	if path.numimages()<3:
		print "A path with less than three images cannot be relaxed.\nUse pathprepare to interpolate some intermediate images!\nAborting."
		sys.exit(1)


	#if forces and energies were successfully read  from checkpoint,
	#previous run is resumed. report that
	stepoffset=path.nstep
	if stepoffset!=0 and options.verbosity>=pypath.constants.VBL_NORMAL:
		print "Loaded checkpoint is a continuation of previous path search after %d steps" % stepoffset
	if path.has_energies() and options.verbosity>=pypath.constants.VBL_NORMAL:
		print "Loaded checkpoint provides initial Energies"
	if path.has_realforces() and options.verbosity>=pypath.constants.VBL_NORMAL:
		print "Loaded checkpoint provides initial Forces"


	startclock=time.clock()
	startwalltime=time.time()
	
	if options.pathdebug:
		if options.verbosity >= pypath.constants.VBL_NORMAL:
			print "Writing path debug output"
		path.writePathDebugInfo=True

	path.nebiterate(scheduler)

	endclock=time.clock()
	endwalltime=time.time()
	
	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print "shutting down scheduler and calculator(s)"
	scheduler.shutdown()

	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print "\nFinished after %d NEB iterations and %d total Energy calculations in %d SCF cycles" % (path.nstep-stepoffset,scheduler.jobCount,scheduler.iterationCount)
		print "Script runtime %f seconds CPU, %f seconds wall time" % (endclock-startclock,endwalltime-startwalltime)
		print "Calculator runtime %f seconds CPU, %f seconds wall time\n\n" % (scheduler.cpuTime,scheduler.wallTime)
##		print "Calculation program runtime %f seconds CPU, %f seconds wall time\n\n" % (calc.workercputime, calc.workerwalltime)

	#some stuff to do when finished
	ratheryoudidnts=["""1. I'd Really Rather You Didn't Act Like A Sanctimonious, Holier-Than-Thou Ass When 
Describing My Noodly Goodness. If Some People Don't Believe in Me, That's Okay. Really, 
I'm Not That Vain. Besides, This Isn't About Them So Don't Change The Subject.""",
	"""2. I'd Really Rather You Didn't Use My Existence As A Means To Oppress, Subjugate, Punish, 
Eviscerate, And/Or, You Know, Be Mean To Others. I Don't Require Sacrifices And Purity Is 
For Drinking Water, Not People""",
	"""3. I'd Really Rather You Didn't Judge People For The Way They Look, Or How They Dress, 
Or The Way They Talk, Or, Well, Just Play Nice, Okay? 
Oh, And Get This Into Your Thick Heads: Woman=Person. Man=Person. Samey-Samey. One Is Not 
Better Than The Other, Unless We're Talking About Fashion And I'm Sorry, But I Gave That To 
Women And Some Guys Who Know The Difference Between Teal And Fuchsia.""",
	"""4. I'd Really Rather You Didn't Indulge In Conduct That Offends Yourself, Or Your Willing,
Consenting Partner Of Legal Age AND Mental Maturity. As For Anyone Who Might Object, I 
Think The Expression Is Fo F*** Yourself, Unless They Find That Offensive In Which Case 
They Can Turn Off The TV For Once And Go For A Walk For A Change.""",
	"""5. I'd Really Rather You Didn't Challenge The Bigoted, Misogynist, Hateful Ideas Of Others 
On An Empty Stomach. Eat, Them Go After The B*******.""",
	"""6. I'd Really Rather You Didn't Build Multimillion-Dollar Churches/Temples/Mosques/Shrines 
To My Noodly Goodness When The Money Could Be Better Spent (Take Your Pick):
  A. Ending Poverty
  B. Curing Deseases
  C. Living In Peace, Loving With Passion And Lowering The Cost Of Cable
I Might Be A Complex-Carbohydrate Omniscient Being, But I Enjoy The Simple Things In Life. 
I Ought To Know. I AM The Creator.""",
	"""7. I'd Really Rather You Didn't Go Around Telling People I Talk To You. You're Not That
Interesting. Get Over Yourself. And I Told You To Love Your Fellow Man, Can't You Take 
A Hint?""",
	"""8. I'd Really Rather You Didn't Do Unto Other As You Would Have Them Do Unto You If You 
Are Into, Um, Stuff That Uses A Lot Of Lubricant/Leather/Las Vegas. If The Other Is Into 
It, However (Pursuant To #4), Then Have At It, Take Pictures, And For The Love Of Mike, 
Wear A CONDOM! Honestly, It's A Piece Of Rubber. If I Didn't Want It To Feel Good When 
You Did IT I Would Have Added Spikes, Or Something."""]
	
	if options.verbosity >= pypath.constants.VBL_NORMAL:
		print random.choice(ratheryoudidnts)
		print "\nFrom Bobby Henderson, The Gospel of the Flying Spaghetti Monster, Villard, New York, 2006\n\n"

	return

def readconfig_options(mode="s"):
	# Default values
	configdefaults=dict(
		maxit="100",
		rmstol="1e-5",
		forcetol="1e-4",
		calculator="dftb",
		fixedatoms="",
		climber="-1",
		fmode='s',
		tangmode='s',
		springk='1.0',
		stepwidth='0.1',
		relmode='s',
		charge='0.0',
	)
	
	# valid calculator types
	# threading mode put on ice for now, as it requires major changes to all calculators
##	VALIDSCHEDS=['s','p','t']
	VALIDSCHEDS=['s','p']

	# Read pypath.ini

	config = ConfigParser.SafeConfigParser(defaults=configdefaults)
	config.read("pypath.ini")
	if not config.has_section("NEB"):
		config.add_section("NEB")
	if not config.has_section("PYPATH"):
		config.add_section("PYPATH")

	#Parse Options

	usage="usage: %prog [options] <startgeometry> [<intermediate geometries>] <endgeometry>"

	parser=OptionParser(usage)

	
	parser.add_option("-m","--maxit",
			action="store", type="int", metavar="N", dest="maxit", default=config.getint("PYPATH","maxit"),
			help="Maximum number of NEB iterations for this run, default=%default")

	parser.add_option("-r","--rmstol",
			action="store", type="float", metavar="RMS", dest="rmstol", default=config.getfloat("PYPATH","rmstol"),
			help="stop iterating if path rms force in a.u. falls below RMS, default=%default")

	parser.add_option("-f","--forcetol", 
			action="store", type="float", metavar="FMAX", dest="forcetol", default=config.getfloat("PYPATH","forcetol"),
			help="stop iterating if path maximum of atomic forces in a.u. falls below FMAX, default=%default")

	parser.add_option("-c","--calculator", 
			action="store", metavar="CALC", dest="calculator", default=config.get("PYPATH","calculator"),
			help="use CALC to calculate energies and forces, default=%default")

	parser.add_option("-x","--fixedatoms", 
			action="store", metavar="LIST", dest="fixedatoms", default=config.get("PYPATH","fixedatoms"),
			help="Keep atoms in whitespace delimeted LIST fixed, default=%default")

	parser.add_option("-s","--stepwidth", 
			action="store", type="float", metavar="D_T", dest="stepwidth", default=config.getfloat("PYPATH","stepwidth"),
			help="Use D_T as stepwidth for velocity Verlet relaxation, default=%default")
	
	parser.add_option("--charge", 
			action="store", type="float", metavar="Q", dest="charge", default=config.getfloat("PYPATH","charge"),
			help="Set Q as the system charge in all path images, default=%default")

	parser.add_option("-v", 
			action="count", dest="verbosity", default=pypath.constants.VBL_NORMAL,
			help="Increase verbosity level, default=%default")
	
	parser.add_option("-q", "--quiet",
			action="store_const", const=pypath.constants.VBL_QUIET, dest="verbosity",
			help="Limit output to fatal errors and critical warnings, no status output at all")
			
	parser.add_option("--silence",
			action="store_const", const=pypath.constants.VBL_SILENCE, dest="verbosity",
			help="Limit output to fatal errors only. Not even warnings.")			

	parser.add_option("-d","--scheduler",
			action="store", type="string", dest="scheduler", default="s",
			help="Choose scheduler to manage execution of single calculations. Choices are 's': serial, 'p': MPI parallel. default=%default")

	parser.add_option("-t","--maxworkthreads",
			action="store", type="int", dest="maxthreads", default=0,
			help="Set maximum number of concurrent parallel threads in thread scheduling mode. 0:= no limit. default=%default")
	
	parser.add_option("--pathdebug",
			action="store_true", default=False, dest="pathdebug",
			help="Store advanced Geometry and forces debug information on path in every iteration")

	# special options for parallel mode go here
	
##	if mode=="p":
       	parser.add_option("-p","",
		action="store", dest="dummy", default=None,
		help="dummy option")


	(options,args) = parser.parse_args()
	

	#special arguments processing for parallel mode goes here
##	if mode=="p":
	# remove some stupid mpich arguments...
	if options.dummy!=None:
		p4wd=args.pop()
		p4pg=args.pop()
	
	#check if at least one mobile image is specified
	if (len(args)<3) and not (len(args)==1 and (args[0].endswith(".fmg") or args[0].endswith(".fmg.gz") or args[0].endswith(".fmg.bz2"))):
		print "You must specify at least three structures. Abort."
		return

	#check if valid calculator type was specified
	if not options.calculator in pypath.constants.PASTACALCS:
		print "Unknown (or unsupported) calculator specified. abort."
		sys.exit(1)
	
	#check if valid scheduler type was specified
	if not options.scheduler in VALIDSCHEDS:
		print "Unknown (or unsupported) scheduler specified. abort."
	
	return (options, config, args)



if __name__ == '__main__':
	# detect MPI environment. If running under MPI, only initialize mpiScheduler as slave
	if os.environ.has_key("LAMRANK"):
		if os.environ["LAMRANK"]=='0':
			mpi=True
			master=True
		else:
			mpi=True
			master=False
	elif os.environ.has_key("MPIRANK"):
		if os.environ["MPIRANK"]=='0':
			mpi=True
			master=True
		else:
			mpi=True
			master=False
	#arminius mpich detection
	elif "-p4amslave" in sys.argv:
		mpi=True
		master=False
	elif "-p4pg" in sys.argv:
		mpi=True
		master=True
	else:
		mpi=False
		master=True
	#master mode
	if master:
		# Store the current path for the signal handlers
		global startpath
		startpath=os.path.realpath(".")
		# Set the signal handlers
		signal.signal(signal.SIGINT, stophandler)
		signal.signal(signal.SIGXCPU, stophandler)
	
		# Parse config file
		nebrun(mode="s")
	elif mpi:
		mpisched=pypath.Schedulers.mpiScheduler(calcinitializer,None,pypath.constants.VBL_NORMAL)
