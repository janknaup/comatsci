#!/usr/bin/env python
##############################################################################
# geoconv
# Part of COmputational MAterials SCIence toolkit - comatsci
# (c) 2005-2008 by Jan M. Knaup <janknaup@gmail.com>
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

def mainfunc():
	import os,sys
	
	from comatsci import geometry,constants
	
	
	from optparse import OptionParser
	
	import random
	
	import numpy
	
	programname=os.path.basename(sys.argv[0]).lower()
	
	geo=geometry.Geometry()
	
	#map mode names, file name endings and writer functions
	workmodes={"xyz"	:	("xyz",geometry.Geometry.writexyz),
		"pdb" :	("pdb",geometry.Geometry.writepdb),
		"fdf" :	("fdf",geometry.Geometry.writefdf),
		"gen" :	("gen",geometry.Geometry.writegen),
		"fmg" : ("fmg",geometry.Geometry.writefmg),
		"xyzq": ("xyzq",geometry.Geometry.writexyzq),
		"tm"  : ("coord",geometry.Geometry.writeTurboMole),
		"aims":	("geometry.in",geometry.Geometry.writeAIMS),
		"cdh" : ("cdh",geometry.Geometry.writeCDH),
		"vasp":	("POSCAR",geometry.Geometry.writecar)}
	
	# geometry formats beloging to the work modes listed here have fixed file names by default.
	# default file name will not be generated from input file name, but the fixed name will be used
	fixednames=("aims","tm","vasp")
	
	# dirty trick to make optparse display an intelligible error message, if output format is not specified
	if not programname[2:] in workmodes.keys():
		programname="togen"
	
	usage="usage: %prog [options] <input file> [<output file>]\nLicensed under the Non-Profit Open Software License version 3.0"
	
	parser=OptionParser(usage)
	
	parser.add_option("-f","--format",
			action="store", metavar="F", default=programname[2:], dest="mode", choices=workmodes.keys(),
			help="Write output geometry in format F. Default=%s. Choose from : %s" % (programname[2:],str(workmodes.keys()),))
	parser.add_option("-x","--extend",
			action="store", metavar="XTND", dest="xtnd", type="string", default="1:1:1",
			help="periodically extend the input geometry to include a:b:c original supercells in the a:b:c lattice directions, default = %default (only the original supercell)")
	##parser.add_option("-c","--coordinate-mode",
	##		action="store", metavar="M", dest="mode", type="string", default="C",
	##		help="write .gen file coordinades in _c_arthesian or _f_ractional mode, default=%default")
	parser.add_option("-p","--population",
			action="store", metavar="F", dest="popfile", type="string", default=None,
			help="Read atomic populations from file F. No default.")
	parser.add_option("-l","--layer",
			action="store", metavar="L", dest="layer", type="int", default=None,
			help="Write only atoms from layer L into output geometry. default: write whole geometry")
	parser.add_option("-e","--element",
			action="store", metavar="E", dest="element", type="int", default=None,
			help="Write only atoms of lement E into output geometry. default: write all elements")
	parser.add_option("-s", "--scale",
					action="store", metavar="S", dest="scale", type="float", default=None,
					help="scale the whole geometry by a factor of S.")
	parser.add_option("--anisoscale",
					action="store", metavar="S", dest="anisoscale", type="float", nargs=3, default=None,
					help="anisotropically scale the whole geometry by a vector S. Scale in lattice vector directions for supercells and Cartesian directions for cluster geometries")
	#.................................................
	parser.set_defaults(append=None)
	parser.add_option("-a","--append",
			action="append", metavar="F L", dest="append", type="string", nargs=2,
			help="Append geometry from file F to input geometry into layer named L. L _must_ be specified, if L ist not present in primary input geometry, it will be created in the output geometry. The special layer name 'default' specifies the .fmg default layer 0.")	
	#.................................................
	parser.set_defaults(lattice=None)
	parser.add_option("-L","--lattice",
			 action="store", dest="lattice", nargs=9, type="float",
			 help="Set mode to periodic and apply supplied lattice vectors. Suppy in format a1 a2 a3 b1 ... c3.")
	parser.set_defaults(foldback=False)
	parser.add_option("-F","--foldback",
			 action="store_true", dest="foldback",
			 help="Fold periodic geometry back into unit cell.")
	parser.add_option("--chargeconstraints",
			action="store_true", dest="chargeconstraints",
			help="Print charge constraints string for specified atoms")
	
	parser.add_option("--tolayer",
			action="store", metavar="L", dest="mtl", type="int", default=None,
			help="Move atoms defined by --atomindices or --atomlist to layer L. L must be defined in the geometry, either --atomindices or --atomlist must be given. Cannot be combined with -l.")
	parser.add_option("--translate",
			action="store", metavar="V", dest="trn", type="string", default=None,
			help="Translate atoms defined by --atomindices or --atomlist by vector V. V must be in the form x:y:z. Either --atomindices or --atomlist must be given.")
	parser.set_defaults(randomvacancies=None)
	parser.add_option("--random-vacancies",
					action="store", metavar="E,N", dest="randomvacancies", type="string",
					nargs=2,
					help="create N vacancies of element E by randomly removing atoms." )
	parser.add_option("--atomindices",
			action="store", metavar="IDX", dest="atindices", type="string", default=None,
			help="Declare, which single atoms are to be modified, argument is a whitespace delimited list of integer atom index numbers. Counts from 0. Mutually exlusive with --atomlist.")
	parser.add_option("--atomlist",
			action="store", metavar="LST", dest="atlist", type="string", default=None,
			help="Declare, which single atoms are to be modified, argument is a whitespace delimited list of integer atom serial numbers. Counts from 1. Mutually exlusive with --atomindices.")
	#.................................................	
	parser.set_defaults(addlayer=None)
	parser.add_option("--addlayer",
			action="store", metavar="L", dest="addlayer", type="string", 
			help="Add layer with name L to geometry.")
	
	(options,args) = parser.parse_args()
	
	
	if (len(args) > 2):
		print "unexpected arguments in command line : %s" % (str(args[2:]))
		print "aborting"
		sys.exit(1)
		
	if (len(args) < 1):
		print "need at least an input file"
		parser.print_usage()
		sys.exit(1)
	
	infilename=os.path.abspath(args[0])
	
	geo.readfile(infilename)
	
	# if a layer to add was specified, do it
	if options.addlayer!=None:
		newLayerIndex=geo.addlayer(options.addlayer)
		print "Layer name '%s' with index %d added to geometry."  %(options.addlayer, newLayerIndex)
	
	# if geometry files to append are specified:
	if options.append !=None:
		# iterate through appended geometry apecifiers
		for i in options.append:
			# be verbose
			print "appending geometry from %s into layer %s ..." %(i[0],i[1]),
			# construct epmty geometry object to append and read the input file
			appendgeo=geometry.Geometry()
			# appended geometry is stored as a tuple (filename,layername)
			appendgeo.readfile(i[0])
			# get layer index from layername
			appendlayer=geo.layerbyname(i[1])
			# if layer does not exist and layer name is 'default', set appendlayer to 0
			if appendlayer==None and i[1]=="default":
				appendlayer=0
			# if layer does not yet exist, create it
			if appendlayer==None:
				appendlayer=geo.addlayer(i[1])
				
			# append geometry
			geo.appendgeometryflat(appendgeo,appendlayer)
			# be verbose
			print "done."
	
	if options.popfile!=None:
		print "reading dftb mulliken populations from %s" % (options.popfile)
		if options.popfile.find("CHR.DAT")!=-1:
			geo.read_dftb_charges(options.popfile)
		elif options.popfile.find("detailed.out")!=-1:
			geo.read_noodle_charges(options.popfile)
		else:
			print "Cannot identify population file type.\nKnown types are:\n\tCHR.DAT\n\tdetailed.out\nabort."
			sys.exit(-1)
	
	if (options.xtnd != "1:1:1"):
		extensions=options.xtnd.split(':')
		if len(extensions)!=3:
			print "periodic extension needs three integer arguments in the form a:b:c. aborting."
			sys.exit(1)
		for i in range(len(extensions)):
			extensions[i]=int(extensions[i])
		geo.periodicexpand(extensions)
	
	if options.layer!=None:
		print "writing layer subgeometry for layer %d:%s" % (
			options.layer,geo.LayerDict[options.layer].Name)
		outgeo=geo.layersubgeometry(options.layer)
	else:
		outgeo=geo
	
	if options.element!=None:
		print "filtering for element %d%s" % (
			options.element,geo.PTE[options.element])
		outgeo=outgeo.elementsubgeometry(options.element)
	else:
		pass
	
	if options.atindices!=None:
		#--atomindices and --atomlist are mutually exclusive!
		if options.atlist!=None:
			print "--atomindices and --atomlist are mutually exclusive! abort."
			sys.exit(1)
		else:
			modatoms=[int(s) for s in options.atindices.split()]
	elif options.atlist!=None:
		modatoms=[ int(s)-1 for s in options.atlist.split()]
		
	else:
		modatoms=None
	
	if modatoms!=None and max(modatoms) >= geo.Atomcount:
		print "nonexistet atom(s) in list of atoms to be modified. abort."
		sys.exit(1)
	
	
	if options.mtl!=None:
		if modatoms==None:
			print "Must specify atoms which should be modified!"
			sys.exit(1)
		elif not outgeo.LayerDict.has_key(options.mtl):
			print "Layer %d, specified in --tolayer does not exist. abort." % (options.mtl)
			sys.exit(1)
		else:
			for i in modatoms:
				outgeo.AtomLayers[i]=options.mtl
	
	
	if options.trn!=None:
		if modatoms==None:
			print "Must specify atoms which should be translated"
			sys.exit(1)
		else:
			tv=[ float(s)/geometry.Angstrom for s in options.trn.split(":")]
			outgeo.translate(tv,modatoms)
	
	if options.chargeconstraints:
		if modatoms==None:
			print "Must specify atoms should be charge-constrained"
			sys.exit(1)
		else:
			print outgeo.getHSDChargeConstraintString(modatoms)
	
	if options.scale!=None:
			print "scaling output geometry by %f " % options.scale
			outgeo.Geometry*=options.scale
			if outgeo.Mode=="S":
					outgeo.Lattice*=options.scale
					
	if options.anisoscale!=None:
		if outgeo.Mode=="S":
			print "anisotropically scaling output geometry along lattice vectors"
			newlattice=numpy.multiply(numpy.array(outgeo.Lattice),numpy.array(options.anisoscale))
			outgeo.Geometry=numpy.dot(outgeo.fractionalGeometry,newlattice)
			outgeo.Lattice=newlattice
		else:
			print "anisotropically scaling output geometry in Cartesian directions"
			outgeo.Geometry=numpy.multiply(outgeo.Geometry,numpy.array(options.anisoscale))
					
	
	if options.lattice!=None:
		print "Applying periodic boundary conditions"
		outgeo.Mode="S"
		outgeo.Lattice=numpy.array(options.lattice).reshape((3,3))/constants.ANGSTROM

	if options.foldback:
		if outgeo.Mode!="S":
			print "Ignoring fold-back request as geometry is not periodic"
		else:
			print "folding coordinates back into unit cell...",
			outgeo.foldToCell()
			print "done."
	
	
	if options.randomvacancies!=None:
		print "Creating random Vacancies:"
		if options.randomvacancies[0].lower() in outgeo.RPTE:
			vacElement=outgeo.RPTE[options.randomvacancies[0].lower()]
		else:
			raise ValueError("Specified Element is unknown: %s"%options.randomvacancies[0])
		vacCount=int(options.randomvacancies[1])
		if outgeo.elementsubgeometry(vacElement).Atomcount < vacCount:
			raise ValueError("More vacancies requested than atoms of element present in geometry.")
		for i in range(vacCount):
			atomSelected=False
			while not atomSelected:
				trialindex=random.choice(range(outgeo.Atomcount))
				if outgeo.AtomTypes[trialindex]==vacElement: atomSelected=True
			outgeo.delatom(trialindex)
		
	if (len(args)==2):
		outfilename=os.path.abspath(args[1])
	elif options.mode in fixednames:
		outfilename=workmodes[options.mode][0]
	else:
		outfilename=os.path.basename(infilename)
		dotpos=outfilename.rfind('.')
		if dotpos==-1:
			dotpos=len(outfilename)
		outfilename=os.path.basename(outfilename[:dotpos]+'.'+workmodes[options.mode][0])
			
	workmodes[options.mode][1](outgeo,outfilename)

if __name__ == '__main__':
	mainfunc()
