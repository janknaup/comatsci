## Automatically adapted for numpy.oldnumeric Oct 27, 2008 by 

#!/usr/bin/python -O

import numpy.oldnumeric as num

A=[-200.,-100.,-170.,15.]
a=[-1.,-1.,-6.5,0.7]
b=[0.,0.,11.,0.6]
c=[-10.,-10.,-6.5,0.7]
x0=[1.,0.,-0.5,-1.]
y0=[0.,0.5,1.5,1.]

xes=num.arrayrange(-2.,1.5,0.05)
yes=num.arrayrange(-1.,2.5,0.05)

outenergies=num.zeros((len(xes),len(yes)),num.Float)

for x in range(len(xes)):
	for y in range(len(yes)):
		# iterate through Mueller-Brown potential terms
		for j in (0,1,2,3):
			#save the exponential value, as we need it for the forces
			expoterm=A[j]*num.exp(a[j]*(xes[x]-x0[j])**2 + b[j]*(xes[x]-x0[j])*(yes[y]-y0[j])+c[j]*(yes[y]-y0[j])**2)
			# store energy contribution
			outenergies[x][y]+=expoterm
			# calculate and store froce contributions
#			tempforces[i][0]+=expoterm*(2.*a[j]*(Geometry.Geometry[i][0]-x0[j])+b[j])
#			tempforces[i][1]+=expoterm*(1.*c[j]*(Geometry.Geometry[i][1]-y0[j])+1.)
		print "%f %f %f" % (xes[x],yes[y],outenergies[x][y])
#	print ""
	print ""
	# finish iterations


