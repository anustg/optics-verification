import numpy as N
from tracer.models.four_parameters_cavity import *
import pickle
from tracer.CoIn_rendering.rendering import *

# -------------------- NRAYS Sensibility analysis ------------------------
#-------------------------Simulation parameters--------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Receiver number of elements (for emissive losses) -------------------------
bins_frustum = 1
bins_cone = 1
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.06 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter
CSR = [0.2]

# An analysis of the precision of the raytrace in function of the number of rays fired.
rayz = 1000
pop = 50000

# Geometries
geoms = N.zeros((5,4))
geoms[0,:] = N.asarray((0.3,0.01,0.05,0.01))
geoms[1,:] = N.asarray((0.3,0.01,0.05,2.99))
geoms[2,:] = N.asarray((0.4,1.5,1.5,1.5))
geoms[3,:] = N.asarray((0.5,0.01,1.5,3.))
geoms[4,:] = N.asarray((0.6, 3., 3., -2.99))

VFs = N.zeros((N.shape(geoms)[0],1+bins_frustum+bins_cone,1+bins_frustum+bins_cone))
areas = N.zeros((N.shape(geoms)[0],1+bins_frustum+bins_cone))

for g in xrange(N.shape(geoms)[0]):
	VF_case = Fourparamcav(geoms[g,0], geoms[g,1], geoms[g,2], geoms[g,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
	VFs[g], areas[g] = VF_case.VF_sim(bins_frustum,bins_cone)


for s in xrange(len(CSR)):
	if CSR[s] == None:
		file_name ='/home/charles/Saves/CSR/AUSES/AUSES_sim_3/Nrayz_sens_AUSES_CSR_None'
	else:
		file_name ='/home/charles/Saves/CSR/AUSES/AUSES_sim_3/Nrayz_sens_AUSES_CSR_%f' %CSR[s]
	data={}

	for g in xrange(N.shape(geoms)[0]):
		# Initialize solution container
		sim_nrays = N.zeros((pop,12))
		for i in xrange(pop):
			case_sens = Fourparamcav(geoms[g,0], geoms[g,1], geoms[g,2], geoms[g,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
			case_sens.VF = VFs[g]
			case_sens.areas = areas[g]
			sim_nrays[i]=case_sens.sim(Tamb=tempAmbiant, Trec=tempReceiver, CSR = CSR[s], nrays=rayz)
		data['geom_%d' %g] = sim_nrays

	data['raytrace_rays']= rayz
	data['total_pop']=pop
	data['VFs']=VFs
	data['areas']=areas
	pickle.dump(data,open(file_name,'w'))



