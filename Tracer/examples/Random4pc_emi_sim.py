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


#-------------------------Simulation parameters--------------------
# Receiver: optics ------------------------------------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Receiver number of elements (for emissive losses) -------------------------
n = 15
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.05 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter


geoms = N.zeros((1,4)) # contains the receiver geometrical parameters
Asm = [] # contains the Assembly object corresponding to the geoms at the same index.

numgeoms = 20000

# Build the geometries array:
index=0 # counter used to dimension the valid geometries array
while index < numgeoms:
	
	apertureRadius = uniform(0.1,1.5) # (m) >0
	apertureDepth = uniform(0.1,2.5) # (m) >abs(coneDepth) & >0
	coneRadius = uniform(0.1,1.5) # (m)	>0
	coneDepth = uniform(-2.5,2.5) # (m) >-apertureRadius
	max_depth = 2.5
	# Check geometries:
	if N.logical_and((apertureDepth+coneDepth) > 0, (apertureDepth+coneDepth) < max_depth):
		geoms.resize((index+1,4))
		geoms[index] = N.asarray((apertureRadius, apertureDepth, coneRadius, coneDepth))

		index=index+1			

# Build the simulation results container for the actual number of rays
emi_sim_results = N.zeros((N.shape(geoms)[0],5))
t_tot = 0
n_error = []
# Simulate all geometries 1 time
for g in range(N.shape(emi_sim_results)[0]):
	t0 = time.clock()
	Asm = Fourparamcav(geoms[g,0], geoms[g,1], geoms[g,2], geoms[g,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
	# Build the Fourparamsim instances with valid geometries.		
	emi_test = Fourparamsim(Asm, tempAmbiant, tempReceiver, n, 100000)
	# Run the simulation for the current number of rays
	emi_sim_results[g] = emi_test.emi_sim()
	t1 = time.clock()-t0
	t_tot = t_tot+t1
	if N.isnan(emi_sim_results[g]).any():
		n_error.append(g)				
	print g+1,' geometries simulated among ',N.shape(geoms)[0],'. ', len(n_error),'VF errors encountered.'
	print 'Simulation time:', t_tot,'s. ', 'Est. remaining time:',t_tot/(60*(g+1))*(N.shape(emi_sim_results)[0]-(g+1)),'min.'

data = {'raytrace_results':emi_sim_results,'n_error':n_error, 'sim_time (s)':t_tot}
pickle.dump(data,open('/home/charles/Saves/emi_sim_4pc','w'))
