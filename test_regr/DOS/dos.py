#!/usr/bin/python
import numpy as num

import comatsci

from scipy.signal import lfilter, correlate

DT=comatsci.DOS.DOS()
#DT.readTaggedOut("results.tag.gz")
DT.readBandOut("band.out.gz")

ORIG=DT.spreadDOS("lorentz",0.01,0.05)
NEW=DT.fastSpreadDOS("lorentz",0.005,0.05)

refdos=DT
refdata=NEW
output=open("newDOS.dat","w")
for s in range(refdos.spins):
    for i in range(len(refdata[0])):
        print >> output, "%12.6f  %12.6f %12.6f %12.6f" % (refdata[0][i],refdata[s+1][0][i],refdata[s+1][1][i],refdata[s+1][2][i])

output.close()


refdata=ORIG
output=open("origDOS.dat","w")
for s in range(refdos.spins):
    for i in range(len(refdata[0])):
        print >> output, "%12.6f  %12.6f %12.6f %12.6f" % (refdata[0][i],refdata[s+1][0][i],refdata[s+1][1][i],refdata[s+1][2][i])

output.close()



