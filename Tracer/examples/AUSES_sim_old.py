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

'''
#-------------------------Simulation parameters---------------------
Simulation file for the AUSES Abstract proposed: Open cavity receiver geometry influence on radiative losses.
'''
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
absDish = 0.06 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter
# Source: sunshape and number of rays ------------------------------
CSR=[None, 0.2, 0.4]
nrays=1000000

numgeoms = 2000
for s in xrange(len(CSR)):
	file_dir = '/home/charles/Saves/AUSES_results_new_%d' %s
	sim_comment= 'absReceiver='+str(absReceiver)+' tempAmbiant='+str(tempAmbiant)+' tempReceiver='+str(tempReceiver)+' n='+str(n)+' dishDiameter='+str(dishDiameter)+' dishFocus='+str(dishFocus)+' absDish'+str(absDish)+' sigma='+str(sigma)+' CSR='+str(CSR[s])+ ' nrays'+str(nrays)
	# Build the simulation results container
	brute_sim_results = N.empty((numgeoms,12))
	t_tot = 0
	g=0
	# Build the geometries array: Random declaration in boundary interval
	while g < numgeoms:
		t0 = time.clock()
		apertureRadius = uniform(0.3,0.6)# (m) >0
		apertureDepth = uniform(0.,3.) # (m) >abs(coneDepth) & >0
		coneRadius = uniform(0.05,3.) # (m) >0
		coneDepth = uniform(-3.,3.) # (m) >-apertureRadius
		max_depth = 3.
		# Check geometries: If the depth is ok, index and go on.
		if N.logical_and((apertureDepth+coneDepth) > 0, (apertureDepth+coneDepth) < max_depth):
			case = Fourparamcav(apertureRadius, apertureDepth, coneRadius, coneDepth, absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
			print apertureRadius, apertureDepth, coneRadius, coneDepth
			brute_sim_results[g] = case.sim(Tamb=tempAmbiant, Trec=tempReceiver, CSR = CSR[s], n=n, nrays=nrays)
			#brute_sim_results[g] = case.ray_sim(CSR = CSR,  nrays=nrays)
			#brute_sim_results[g] = case.emi_sim(Tamb=tempAmbiant, Trec=tempReceiver, n=n)
			g+=1
	
			t_tot = t_tot+time.clock()-t0

		if g%20. == 0.:
			save = open(file_dir,'w')
			pickle.dump({'raytrace_results':brute_sim_results, 'sim_time (s)':t_tot, 'sim_comment':sim_comment},save)	
			save.close()
				
		print g+1,' geometries simulated among ',numgeoms,'.'

		print 'Simulation time:', t_tot,'s. ', 'Est. remaining time:',t_tot/(60*(g+1))*(N.shape(brute_sim_results)[0]-(g+1)),'min.'

	data = {'raytrace_results':brute_sim_results, 'sim_time (s)':t_tot, 'sim_comment':sim_comment, 'raytrace_rays':nrays}
	save_end = open(file_dir,'w')
	pickle.dump(data,save_end)
	save_end.close()
	print file_dir
	print sim_comment
