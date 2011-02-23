#!/usr/bin/python

import comatsci

geo=comatsci.Geometry.AnalysisGeometry()
geo.readfile("alpha_quartz_ideal_1x1x1.gen")

geo.locateVacancies(specvalences={8:2,14:4})

geo=comatsci.Geometry.AnalysisGeometry()
geo.readfile("alpha_quartz_1Vac_1x1x1.gen")

vg=geo.locateVacancies(specvalences={8:2,14:4})

print vg.xyzstring()

geo=comatsci.Geometry.AnalysisGeometry()
geo.readfile("start.gen")

vg=geo.locateVacancies(specvalences={8:3,22:6})

print vg.xyzstring()

geo=comatsci.Geometry.AnalysisGeometry()
geo.readfile("geo_end.gen.gz")

vg=geo.locateVacancies(specvalences={8:3,22:6},tolerance=1.3)

print vg.xyzstring()


