import numpy as N
import matplotlib.pyplot as plt
import pickle
from random import uniform

from tracer.models.four_parameters_cavity import *
from tracer.surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.tracer_engine import *
from tracer.CoIn_rendering.rendering import *

#-------------------------Simulation parameters--------------------
# Receiver: optics ------------------------------------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Receiver number of elements (for emissive losses) -------------------------
n = 10
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.05 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter

geoms = N.zeros((1,4)) # contains the receiver geometrical parameters

numgeoms = 4000

# Build the geometries array:
index=0 # counter used to dimension the valid geometries array
while index < numgeoms:
	
	apertureRadius = uniform(0.3,0.5) # (m) >0
	apertureDepth = uniform(0.05,3.) # (m) >abs(coneDepth) & >0
	coneRadius = uniform(0.05,3.) # (m)	>0
	coneDepth = uniform(-3.,3.) # (m) >-apertureRadius
	max_depth = 3.
	# Check geometries:
	if N.logical_and((apertureDepth+coneDepth) > 0, (apertureDepth+coneDepth) < max_depth):
		geoms.resize((index+1,4))
		geoms[index] = N.asarray((apertureRadius, apertureDepth, coneRadius, coneDepth))

		index=index+1					

# Build the simulation results container
brute_sim_results = N.zeros((N.shape(geoms)[0],12))
t_tot = 0
n_error = []
# Simulate all geometries 1 time
for g in range(N.shape(brute_sim_results)[0]):
	t0 = time.clock()
	case = Fourparamcav(geoms[g,0], geoms[g,1], geoms[g,2], geoms[g,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
	# Run the simulation
	brute_sim_results[g] = case.sim(Tamb=tempAmbiant, Trec=tempReceiver, CSR = 0.2, n=n, nrays=1000000)
	del(case)
	t1 = time.clock()-t0
	t_tot = t_tot+t1
	if g%10 == 0:
		data = {'raytrace_results':brute_sim_results,'n_error':n_error, 'sim_time (s)':t_tot}
		pickle.dump(data,open('/home/charles/Saves/Random_brute_sim_1000000_15_CSR02','w'))	
	if N.isnan(brute_sim_results[g]).any():
		n_error.append(g)				
	print g+1,' geometries simulated among ',N.shape(geoms)[0],'. ', len(n_error),'VF errors encountered.'
	print 'Simulation time:', t_tot,'s. ', 'Est. remaining time:',t_tot/(60*(g+1))*(N.shape(brute_sim_results)[0]-(g+1)),'min.'
data = {'raytrace_results':brute_sim_results,'n_error':n_error, 'sim_time (s)':t_tot}
pickle.dump(data,open('/home/charles/Saves/Random_brute_sim_1000000_15_CSR02','w'))
print 'File: /home/charles/Saves/Random_brute_sim_1000000_15_CSR02'
