import numpy as N
import matplotlib.pyplot as plt
import pickle

from tracer.models.four_parameters_cavity import *
from tracer.surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.tracer_engine import *


# timer start
t0 = time.clock()

#-------------------------Simulation parameters--------------------
# Receiver: 4 parameters cavity geometry ---------------------------------
apertureRadius = N.arange(0.2,1.,0.2) # (m) >0
apertureDepth = N.arange(0.2,2.,0.2) # (m) >abs(coneDepth)
coneRadius = N.arange(0.2,1.,0.2) # (m)	>0
coneDepth = N.arange(-1.,1.,0.2) # (m) >-apertureRadius
# Receiver: optics ------------------------------------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 650+273.15 # Receiver internal wall temperature in K
# Receiver number of elements (for emissive losses) -------------------------
n = 20
# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.4 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.05 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter

'''
_________________________________________________________________________________________________
Simulation optimisation routines:

class Fourparamcav():
	def __init__(self, apertureRadius, apertureDepth, coneRadius, coneDepth, absReceiver, emsReceiver,dishDiameter, dishFocus, absDish, sigma):

class Fourparamsim():
	def __init__(self, Asm, tempAmbiant, tempReceiver, n, nrays):

	Methods called: 
	- self.ray_sim(): raytrace of reflective losses only.
		returns: self.apertureRadius, self.apertureDepth, self.coneRadius, self.coneDepth, self.flux_input, self.Ref_losses_tot, self.ReceiverA
	- self.emi_sim(): emissive losses evaluation only.
		returns: self.emissive_losses (total emissive losses through the aperture).
	- self.sim(): combined simulation of emissive and reflective losses
		returns: self.apertureRadius, self.apertureDepth, self.coneRadius, self.coneDepth, self.total_losses, self.emi_losses, self.RealReceiverA, self.System_efficiency

_________________________________________________________________________________________________
'''

'''

# Simple cases...................................................................................................................................................................................................................................................

cavity1 = Fourparamsim(0.5, 0.9 , 0.6, 0.1, absReceiver, emsReceiver, tempAmbiant, tempReceiver, n, dishDiameter, dishFocus, absDish, sigma, 1000000)
cavity2 = Fourparamsim(0.5, 0.9 , 0.7, 0.1, absReceiver, emsReceiver, tempAmbiant, tempReceiver, n, dishDiameter, dishFocus, absDish, sigma, 1000000)

cavity1_results = cavity1.sim()
cavity2_results = cavity2.sim()

fluxplot=plt.figure()
for test in [cavity1,cavity2]:
	# fluxmap for a specific case:
	if test == cavity1:
		plt.subplot(121)
	else:
		plt.subplot(122)
	test.fluxmap()	

plt.show()
'''

# Combined raytrace and geometrical optimization:
# 1- We sim() 1 time each case at nrays=1000
# 2- Take the minimum and run it until its losses average estimator's standard deviation gets below a determined threshold (5% of the average estimator?). Compute the standard deviation of the losses.
# 3- Test all points being under the average estimator  of the range and store their geometries.
# 4- Take the non-stored geometries and compute a representative spread of results on them 
#	  to ensure that they do not belong to the minimum point spread. If they do, take their
# 	  geometry and add it to the store. To compute this spread, the tracer runs eac
# 5- If there is a new global minimum (lower than the lower value of the minimum spread),
#    take this as the new minimum and start again from 2-.
# 6- Set nrays to the next value and start again from 1 with the new set of geometries.
 
spread = 5 # Number of times the minimum is calculated to get the dispersion of results of this geometry.
geoms = N.zeros((1,4)) # contains the receiver geometrical parameters
bad_geoms = []
Asm = [] # contains the Assembly object corresponding to the geoms at the same index.
rayz = N.array([500, 1000, 5000, 10000, 50000])

# Build the geometries array:
index=0 # counter used to dimension the valid geometries array
for ar in range(len(apertureRadius)):
	for ad in range(len(apertureDepth)):
		for cr in range(len(coneRadius)):
			for cd in range(len(coneDepth)):
				# Check geometries:
				if (apertureDepth[ad]+coneDepth[cd]) > 0:
					geoms.resize((index+1,4))
					geoms[index] = N.asarray((apertureRadius[ar], apertureDepth[ad], coneRadius[cr], coneDepth[cd]))
					Asm.append(Fourparamcav(apertureRadius[ar], apertureDepth[ad], coneRadius[cr], coneDepth[cd], absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma))
					
					index=index+1					
				else:	
					bad_geoms.append((apertureRadius[ar], apertureDepth[ad], coneRadius[cr], coneDepth[cd]))
					continue

print 'Valid geoms:', N.shape(geoms)[0]

for r in range(len(rayz)):
	print 'Number of rays: ', rayz[r] 

	# Build the simulation results container for the actual number of rays
	brute_sim_results = N.zeros((N.shape(geoms)[0],8))
	
	# Simulate all geometries 1 time
	for g in range(N.shape(brute_sim_results)[0]):
		# Build the Fourparamsim instances with valid geometries.		
		brute_test = Fourparamsim(Asm[g], tempAmbiant, tempReceiver, n, rayz[r])
		# Run the simulation for the current number of rays
		brute_sim_results[g] = brute_test.sim()						
		print g+1,' geometries simulated among ',N.shape(geoms)[0], '. \r'
		assert ~N.isnan(brute_sim_results[g]).any()
			
	
	# Identify the minimum total losses:
	brute_min = N.argmin(brute_sim_results[:,4])
	print 'brute_min: ', brute_min

	# Initialise solution container array for min case.
	brute_min_case_results = N.zeros((spread,7))
	# Compute spread range of the minimum total losses case:
	brute_min_case = Fourparamsim(Asm[brute_min], tempAmbiant, tempReceiver, n, rayz[r])
	while 
		brute_min_case_results[i] = brute_min_case.ray_sim()
		# Calculate total_losses average statistical estimator and its variance
		brute_min_case_avg_est = 1/(i+1)*N.sum(brute_min_case_results[:i,5])
		if i == 0:
			continue
		var = 1/i*N.sum(brute_sim_case[:i,5]**2)-(i+1)/i*brute_min_case_avg_est**2
		std_dev = N.sqrt(var)
		
		

	# Get the max value of this range (take the results from this ray_sim and add the emissive losses calculated previously for this case):
	max_range = N.amax((brute_min_case_results[:,5]))+brute_sim_results[brute_min,5]
	# Store all geometries with total losses under the max_range value:
	efficient = brute_sim_results[:,4] <= max_range
	print 'efficient geometries (in min range): ', N.sum(efficient)
	print 'inefficient geometries (in min range): ', N.sum(~efficient)
	store = N.argwhere(efficient)
	print 'store', N.shape(store)
	print store

	# Take non-stored cases:
	non_store = N.argwhere(~efficient)
	# Initialise solution container array.
	bad_geoms_case_results = N.zeros((spread,7))
	# Loop through non_store cases until they fit in the range or depletion of the spread counter. If they fit, add them to the store and continue the loop; otherwise discard them.
	# TODO, could be done better using tests on standard deviation or an estimator and its variation through the successive raytraces. If the results of the bad_geom_case is expected to be over max_range, it is taken out.
	for ns in range(len(non_store)):
		bad_geoms_case = Fourparamsim(Asm[non_store[ns]], tempAmbiant, tempReceiver, n, rayz[r])		
		for j in range(spread):	
			bad_geoms_case_results[j] = bad_geoms_case.ray_sim()
			# Test total losses result at each iteration:
			if (bad_geoms_case_results[j,5]+brute_sim_results[non_store[ns],5]) <= max_range:
				store.append[non_store[ns]]
				print 'j: ',j 
				print 'store', N.shape(store)
				N.delete(non_store, ns, 0)				
				break						
		print ns+1,' bad_geoms_cases treated over ', len(non_store), '\r'

	# Update geometries to test for next iteration:
	geoms = geoms[store].reshape(len(store), N.shape(geoms)[1])	
	Asm = Asm[store]
	if N.shape(geoms)[0]<=2:
		print 'Asm', N.shape(Asm), Asm		
		print 'geoms', N.shape(geoms)
		print geoms
		break

t1 = time.clock()-t0
print 'total calculation time: ',t1,'s'


'''
# Basic simulations

# Simulation over APERTURE range using constant depth(1m), frustum and cone half angles(60deg)
simaperture = N.empty((len(apertureDiameter),6))
testa = []
for a in range(len(apertureDiameter)):
	testa.append(fourparamsim(apertureDiameter[a], 1, frustumAngle[11], coneAngle[11], absReceiver, emsReceiver, tempAmbiant, tempReceiver, n, dishDiameter, dishFocus, absDish, sigma,1000))
	resulta = testa[a].sim()
	simaperture[a] = N.array(resulta)

# Simulation over DEPTH range using 0.7m aperture diameter, constant frustum angle (60deg) and evolving cone angle.
simdepth = N.empty((len(depth),6))
testd = []
for d in range(len(depth)):
	testd.append(fourparamsim(0.7, depth[d], frustumAngle[11], N.arctan(1.4/depth[d]), absReceiver, emsReceiver, tempAmbiant, tempReceiver, n, dishDiameter, dishFocus, absDish, sigma))
	resultd = testd[d].sim()
	simdepth[d] = N.array(resultd)
'''
'''
______________________________________________________
plots
______________________________________________________
'''
# Brute force low rays number test:
'''
eff = N.arange(len(efficient_sim))
plt.plot(eff, efficient_sim[:,5])


# Efficiency vs. parameters plots:
plots = plt.figure()
ap_plot,dep_plot = plots.add_subplot(2,1,1), plots.add_subplot(2,1,2)
ap_plot.plot(apertureDiameter, simaperture[:,5])
ap_plot.set_title('Efficiency vs. aperture (depth=1m, frustumA=ConeA=60deg)')
dep_plot.plot(depth, simdepth[:,5])
dep_plot.set_title('Efficiency vs. depth (apertureD=0.7m, frustumA=ConeA=60deg)')


plt.show()
'''
