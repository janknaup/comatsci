#!/usr/bin/python

import comatsci
g=comatsci.Geometry.Geometry()
g.readAimsFile("geometry.in")
g.writexyz("aimstest.xyz")
g.writegen("aimstest.gen")
g.writeAIMS("aimstest.in")
