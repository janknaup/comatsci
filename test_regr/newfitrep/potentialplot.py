#!/usr/bin/python
##############################################################################
# potentialplot
# Part of comatsci computational materials science toolkit
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
# output r,U(r),U'(r) for a given pair potential file
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################


#@author: Jan M. Knaup <janknaup@gmail.com>
#@organization: Bremen Center for Compuational Materials Science
#@license: Open Software License version 3.0
#@copyright: Jan M. Knaup  <janknaup@gmail.com>


from comatsci.Calculators import Potentials
import os,sys

helpMessage="""potentialtest (c) Jan M. Knaup 2011
print r,U(r),U'(r) of a given potenitalfunction within the defined cutoff ranges
usage:
potentialtest <potentialfile> <number of points> [<inner cutoff> [outer cutoff]]"""

# try to parse command line options, don?t overburden this, usage should be
# almost foolproof!

if len(sys.argv)<3:
	print >> sys.stderr, "ERROR: missing required arguments."
	print >> sys.stderr, helpMessage
	sys.exit(1)

fileName=sys.argv[1]
if not os.path.exists(fileName):
	print >> sys.stderr, "ERROR: input file does not exist."
	print >> sys.stderr, helpMessage
	sys.exit(1)
else:
	try:
		potential=Potentials.getPotentialFromFile(fileName)
	except:
		print >> sys.stderr, "ERROR: potential could not be read from input file."
		print >> sys.stderr, helpMessage
		raise

try:
	npoints=int(sys.argv[2])
except:
	print >> sys.stderr, "ERROR: number of points could not be understood."
	print >> sys.stderr, helpMessage
	raise

(innercut,outercut)=potential.getCutoffs()


#use semireasonable default values
if innercut==None:
	innercut=0.0000000001
	print >> sys.stderr, "WARNING: setting default inner cutoff %d" % innercut

if outercut==None:
	outercut=10.0
	print >> sys.stderr, "WARNING: setting default outer cutoff %d" % outercut
	

if len(sys.argv)==4:
	try:
		innercut=float(sys.argv[3])
	except:
		print >> sys.stderr, "ERROR: specified inner cutoff could not be understood"
		print >> sys.stderr, helpMessage
		raise

if len(sys.argv)==5:
	try:
		outercut=float(sys.argv[4])
	except:
		print >> sys.stderr, "ERROR: specified outer cutoff could not be understood"
		print >> sys.stderr, helpMessage
		raise

# iterate through ranges and output

step=((outercut-innercut)/float(npoints-1))
for i in range(npoints):
	r=innercut+i*step
	print "%f\t%f\t%f"%(r,potential.value(r),potential.derivative(r))

# done.