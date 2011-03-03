import comatsci
g=comatsci.Geometry.AnalysisGeometry()
r=comatsci.Geometry.AnalysisGeometry()
g.readfile("geo_end.gen.gz")
r.readfile("reference.gen")
v=g.groupVacanciesByReference(r,style="distsort",interstitialLayer=True)
print v.xyzstring()
v.writexyz("vac.xyz")
v.writepdb("vac.pdb")
