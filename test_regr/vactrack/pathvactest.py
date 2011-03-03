import comatsci


#g=comatsci.Geometry.AnalysisGeometry()
#g.readfile("start.gen")
#vg=g.locateVacancies(neighborMethod="neighborTypes",groupMethod="distance",specvalences={22:6,8:3},tolerance=1.2,canonicalNeighbors={22:{8:6},8:{22:3}})

#ng=g.vacanyNeighborsByBondPartners(canonicalNeighbors={22:{8:6},8:{22:3}})


path=comatsci.Path.Reactionpath('checkpoint',[],'s',0.1,0.1,1)

refgeo=comatsci.Geometry.Geometry()
refgeo.readfile("reference.gen")

path.readXyzPath("geo_end.xyz",comatsci.Geometry.AnalysisGeometry)
path.Lattice=refgeo.Lattice
path.Mode=refgeo.Mode


for i in range(path.numimages()):
 path.geos[i].Lattice=refgeo.Lattice
 path.geos[i].Mode=refgeo.Mode

#path.movingAverageResample(50)

path.writexyzpath("smoothed.xyz")

vpath=comatsci.Path.Reactionpath('checkpoint',[],'s',0.1,0.1,1)

refgeo=comatsci.Geometry.Geometry()
refgeo.readfile("reference.gen")

vac=path.geos[0].locateVacancies(neighborMethod="none",groupMethod="reference_b",
		specvalences={22:6,8:3},tolerance=0.49,
		canonicalNeighbors={22:{8:6},8:{22:3}},refold=True,reference=refgeo)
lastvac=vac
initvaccount=int(vac.Atomcount)
vpath.appendGeoObject(vac,checkCompat=False)

for i in range(1,path.numimages()):
	vac=path.geos[i].locateVacancies(neighborMethod="none",groupMethod="reference_b",
		specvalences={22:6,8:3},tolerance=0.49,
		canonicalNeighbors={22:{8:6},8:{22:3}},refold=True,reference=refgeo)

	if vac.Atomcount!=initvaccount:
		vpath.appendGeoObject(lastvac,checkCompat=False)
	else:
		vpath.appendGeoObject(vac,checkCompat=False)
		lastvac=vac

print vac.xyzstring()
vpath.writexyzpath("vac.xyz")


