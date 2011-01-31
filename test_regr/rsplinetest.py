import comatsci
path=comatsci.Path.Reactionpath('checkpoint',[],'s',0.1,0.1,1)
path.readfmgpath("preppath.fmg")
path._genRsplineRep()
parms=[]
for i in range(21):
  parms.append((path._rSplineRep.totalLength/20)*i)

path.rennerSubsplineResample(parms)
