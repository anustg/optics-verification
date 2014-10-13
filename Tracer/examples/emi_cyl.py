import numpy as N
import matplotlib.pyplot as plt
from random import uniform

import pickle

from tracer.models.six_parameters_cylinder import *
from tracer.models.six_parameters_cylinder import *

# Simulation parameters:

apradius = 0.083 #(m)
cylradius = apradius*N.arange(1.1,3.1,0.01)
cyldepth = apradius*N.arange(0.1,2.,0.01)
emissivity = 0.87
n = 1000

tempambiant = 20+273.15
tempreceiver = [450+273.15,550+273.15,650+273.15]

results = N.zeros((5,(len(cyldepth)*len(cylradius)+1)*len(tempreceiver))) # apreture, radius, depth, tempreceiver, losses
index=0
for t in tempreceiver:
	# Flat plate:
	flatplate = FlatplateRec(apradius, emissivity, n)
	losses, tempext, temprec = flatplate.emi_sim(tempambiant, t)
	results[:,index]=N.asarray([apradius,apradius, 0, t, losses])
	index=index+1
	# Same aperture but diameter or depth:
	for r in cylradius:
		for d in cyldepth:
			sixpcyl = Sixparamcyl(apradius, r, d, emissivity, n)
			losses, tempext, temprec = sixpcyl.emi_sim(tempambiant, t)
			results[:,index] = N.asarray([apradius,r, d, t, losses])
			index=index+1

'''
#geoms = N.zeros((len(OR)*len(AR),2))
geoms = N.zeros((len(OR)*ngeoms,2))
index=0
while index<len(OR)*ngeoms:
	gu = uniform(0.005,3.5)*cylradius
	for i in OR:
		geoms[index,1] = gu
		geoms[index,0] = i*cylradius #uniform(0.,1.)*cylradius
		index +=1

for a in OR:
	aperture = a*cylradius
	for d in AR:
		cyldepth = d*cylradius
		geoms[index] = N.asarray([aperture,cyldepth])
		index+=1

print N.shape(geoms)
sixpcyl_results = N.zeros((ngeoms*len(tempreceiver),8))	
for g in xrange(ngeoms):
	cyl = Sixparamcyl(geoms[g,0], cylradius, geoms[g,1], emissivity, n)
	for t in xrange(len(tempreceiver)):
		#print 'geoms:', geoms[g]
		losses, Ta, Tr = cyl.emi_sim(tempambiant, tempreceiver[t])
		sixpcyl_results[g*len(tempreceiver)+t] = N.asarray((geoms[g,0], cylradius, geoms[g,1], emissivity, n, losses, Ta, Tr))
'''
#pickle.dump(results, open('/home/charles/Saves/emi_cyl_sim_apcomp','w'))
