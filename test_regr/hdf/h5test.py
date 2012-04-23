import comatsci,numpy,h5py,cProfile
g=comatsci.Geometry.Geometry()
g.readfile("15qmz.fmg")
g.writeCDH("test.cdh")

path=comatsci.Path.NEBPath('',[],1.,'s','w',-1,1.,'s',1.,1.,1,'v',0,comatsci.constants.VBL_DEBUG1)
path.readXyzPath("geo_end.xyz")
lattice=numpy.array([[  0.6275120000E+01  ,  0.6275120000E+01   , 0.0000000000E+00], [  -0.9187460000E+01   , 0.9187460000E+01 ,   0.0000000000E+00],[    0.0000000000E+00  ,  0.0000000000E+00  ,  0.9900000000E+02]],dtype=float)

for i in range(path.numimages()):
  path.geos[i].Mode="S"
  path.geos[i].Lattice=numpy.array(lattice)

cProfile.run('path.writeCDHPath("testpath.cdh")')

cProfile.run('path.writeCDHPath("testpath-savespace.cdh",savespace=True)')

t=comatsci.Geometry.Geometry()
t.readCDHFile("test.cdh")

tpath=path=comatsci.Path.NEBPath('',[],1.,'s','w',-1,1.,'s',1.,1.,1,'v',0,comatsci.constants.VBL_DEBUG1)
tpath.readCDHPath("testpath.cdh")

tpath2=path=comatsci.Path.NEBPath('',[],1.,'s','w',-1,1.,'s',1.,1.,1,'v',0,comatsci.constants.VBL_DEBUG1)
tpath2.readCDHPath("testpath-savespace.cdh")


