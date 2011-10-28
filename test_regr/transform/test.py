#!/usr/bin/python

import comatsci

g=comatsci.Geometry.Geometry()
g.readfile("h2odimer.xyz")

g.centerOnAtom(4)
g.writexyz("centered.xyz")

anglesteps=20
stepangle=2.0*comatsci.constants.PI*(1.0/anglesteps-1)

g.writexyz("rotated-000.xyz")

for i in range(1,anglesteps):
  print i
  g.rotateAxis("z",stepangle,(3,5))
  filename="rotated-%03i.xyz"%i
  g.writexyz(filename)


