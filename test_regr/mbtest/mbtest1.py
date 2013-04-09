## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

#!/usr/bin/python

import pypath as pp
import numpy.oldnumeric as num

c=pp.Calculators.muellerBrownCalc(verbosity=9)
g=pp.Geometry.Geometry()
#g.addatom(1,(-0.5,1.4,0.))
g.addatom(1,(0.6,0.,0.))
#g.addatom(1,(0.0,-0.05,0.467))
c.start(g,"bla")
e=c.getenergy()
f=c.getforces()

import pypath.Optimizers as opt
optim=opt.steepestDescentOptimizer({"verbosity":9,"maxIterations":10000,"maxF":0.001,"stepSize":0.00001})
#optim=opt.newtonRaphsonOptimizer({"verbosity":9,"maxIterations":10000,"maxF":0.001,"stepSize":0.00001})
#optim=opt.velocityVerletOptimizer({"verbosity":9,"maxIterations":10000,"maxF":0.1,"stepSize":0.00001,"projectVelocities":True})
#optim=opt.velocityVerletOptimizer({"verbosity":9,"maxIterations":10000,"maxF":0.001,"stepSize":0.00001,"projectVelocities":True,"adaptive":True})


c.start(g,"bla")
e=c.getenergy()
f=c.getforces()
print e
print f
g.writexyz("traj.xyz")
while not optim.converged:
	print "step %d" % optim.iterations
	c.start(g,"bla")
	e=c.getenergy()
	f=-c.getforces()
	gshape=num.shape(g.Geometry)
	newx=optim.optStep(g.Geometry.ravel(),e,f.ravel())
	g.Geometry=num.reshape(newx,gshape)
	g.writexyz("traj.xyz","a")

print e
print f
print g.Geometry

