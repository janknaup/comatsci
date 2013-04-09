#!/usr/bin/python
import pypath
import sys

path=pypath.Path.Reactionpath("",[],"s",0.1,0.1,1,0,pypath.constants.VBL_SILENCE)
path.readfmgpath(sys.argv[1])

for i in range(path.numimages()):
	if path.has_energies():
		print "%f\t%f\t%f\t%f" % (path.geos[i].Geometry[0][0],path.geos[i].Geometry[0][1],path.geos[i].Geometry[0][2],path.energies[i])
	else:
		print "%f\t%f\t%f" % tuple(path.geos[i].Geometry[0])

print "\n"
