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
bins_frustum = 1
bins_cone = 1
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.06 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter
# Source: sunshape and number of rays ------------------------------
CSR = 0.2
nrays = 1000000
# Simulation routine parameters:
numgeoms = 1000
# Data from sensibility studies to estimate raytrace accuracy:
rays = [1,2,3,4,5,10,15,20,30,40]
IC = [0.039, 0.028,0.023,0.019,0.017,0.012,0.01,0.009,0.007,0.006]

for r in xrange(len(rays)):
	if r!=0: 
		# Load and filter datas from the previous raytrace according to the simulated IC boundary if this is not the first simulation of the routine.
		file_in = '/home/charles/Saves/CSR/AUSES/AUSES_sim_3/AUSES_results_re_%de6' %rays[r-1]
		data = pickle.load(open(file_in,'r'))
		base_data = data['raytrace_results']
		VFs = data['VFs']
		areas = data['areas']
		good_data = base_data[N.where((base_data[:,5]-base_data[:,6])>(1-IC[r])*N.amax(base_data[:,5]-base_data[:,6]))]
		VFs = VFs[N.where((base_data[:,5]-base_data[:,6])>(1-IC[r])*N.amax(base_data[:,5]-base_data[:,6]))]
		areas = areas[N.where((base_data[:,5]-base_data[:,6])>(1-IC[r])*N.amax(base_data[:,5]-base_data[:,6]))]
		# Actualise number of rays
		rays_so_far = pickle.load(open(file_in,'r'))['raytrace_rays']

	else:
		geoms = N.zeros((1,4)) # contains the receiver geometrical parameters
		# Build the geometries array:
		index = 0 # counter used to dimension the valid geometries array
		while index < numgeoms:
			apertureRadius = uniform(0.3,0.6) # (m) >0
			apertureDepth = uniform(0.,3.) # (m) >abs(coneDepth) & >0
			coneRadius = uniform(0.05,3.) # (m)	>0
			coneDepth = uniform(-3.,3.) # (m) >-apertureRadius
			max_depth = 3.
			# Check geometries:
			if N.logical_and((apertureDepth+coneDepth) > 0, (apertureDepth+coneDepth) < max_depth):
				geoms.resize((index+1,4))
				geoms[index] = N.asarray((apertureRadius, apertureDepth, coneRadius, coneDepth))
				index=index+1			
		good_data = N.zeros((N.shape(geoms)[0],12))
		good_data[:,:4] = geoms
		rays_so_far = 0
		VFs = N.zeros((numgeoms,1+bins_frustum+bins_cone,1+bins_frustum+bins_cone))
		areas = N.zeros((numgeoms,1+bins_frustum+bins_cone))
		for i in xrange(N.shape(good_data)[0]):
			VF_case = Fourparamcav(good_data[i,0], good_data[i,1], good_data[i,2], good_data[i,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
			VFs[i], areas[i] = VF_case.VF_sim(bins_frustum,bins_cone)

	file_out = '/home/charles/Saves/CSR/AUSES/AUSES_sim_3/AUSES_results_re_%de6' %rays[r]

	# Build the new simulation results container
	new_sim_results = N.zeros((N.shape(good_data)[0],12))
	sim_results = N.zeros((N.shape(good_data)[0],12))
	t_tot = 0

	# Calculate the number of 'nrays' bundles necessary to get to rays[r] rays at this stage:
	passes = rays[r]-rays_so_far/nrays

	# Run news sims
	for i in xrange(N.shape(good_data)[0]):
		t0 = time.clock()
		pass_results = N.zeros((passes,12))
		for j in xrange(passes):
			case = Fourparamcav(good_data[i,0], good_data[i,1], good_data[i,2], good_data[i,3], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
			case.VF = VFs[i]
			case.areas = areas[i]
			pass_results[j,:] = case.sim(Tamb=tempAmbiant, Trec=tempReceiver, CSR = CSR, nrays=nrays)
		new_sim_results[i,:4] = good_data[i,:4]
		new_sim_results[i,4:] = N.sum(pass_results[:,4:], axis=0)/passes
		
		# Actual number of rays for the simulation:		
		rays_total = rays_so_far+passes*nrays

		# Timer:
		t_tot = t_tot+time.clock()-t0
	
		# Estimator of the result with the new number of rays
		sim_results[i,:4] = new_sim_results[i,:4]
		sim_results[i,4:] = (good_data[i,4:]*rays_so_far+new_sim_results[i,4:]*passes*nrays)/rays_total
	 	# Comment for the sim dict.
		sim_comment= 'absReceiver='+str(absReceiver)+' tempAmbiant='+str(tempAmbiant)+' tempReceiver='+str(tempReceiver)+' dishDiameter='+str(dishDiameter)+' dishFocus='+str(dishFocus)+' absDish'+str(absDish)+' sigma='+str(sigma)+' CSR='+str(CSR)+ ' nrays'+str(rays_total)
		# Save results al along the raytrace
		if i%20. == 0.:
			save = open(file_out,'w')
			pickle.dump({'raytrace_results':sim_results, 'sim_time (s)':t_tot, 'sim_comment':sim_comment,'raytrace_rays':rays_total, 'VFs':VFs, 'areas':areas}, save)	
			save.close()

		print i+1,' geometries simulated among ',N.shape(good_data)[0],'.'
		print 'Simulation time:', t_tot,'s. '


	# Final save for this number of rays
	data = {'raytrace_results':sim_results, 'sim_time (s)':t_tot, 'sim_comment':sim_comment, 'raytrace_rays':rays_total, 'VFs':VFs, 'areas':areas}
	save_end = open(file_out,'w')
	pickle.dump(data,save_end)
	save_end.close()
	print file_out
	print sim_comment
