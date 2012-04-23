#!/usr/bin/python
import numpy as num

import comatsci

from scipy.signal import lfilter, correlate

DT=comatsci.DOS.DOS()
#DT.readTaggedOut("results.tag.gz")
DT.readBandOut("band.out.gz")

L=DT.lorentzDOS()

J=correlate(L[1][1],L[1][2],"full")

step=L[0][1]-L[0][0]

delta=0.0
for i in range(len(J)/2,0,-1):
 print delta, "\t",J[i]
 delta+=step

