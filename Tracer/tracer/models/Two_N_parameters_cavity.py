'''
__________________________________________________________________________________________________
*********** 2 N Parameters cavity receiver  ************
__________________________________________________________________________________________________
'''
import numpy as N
import time
import matplotlib.pyplot as plt
import matplotlib.cm as cm
N.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.quadric import *
from tracer.paraboloid import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *

from freesteam import *

from view_factors_3D import *
from emissive_losses import *

class TwoNparamcav():
	'''
	TwoNparamcav object initializes the scene asembly for the geometrical and optical parameters given.
	Several methods can be called to perform raytracing and emissive losses calculations.
	- self.ray_sim(): performs a raytrace and returns the effective input flux, the total reflective losses and the energy absorbed by the receiver.
	- self.temperature_guess(): Evaluates the temperature of the tubes considering water/steam flow conditions, incident ratiation and thermal emissions.
	- self.emi_sim(): performs an emissive losses calculation and returns the total emissive losses.
	- self.VF_sim(): Receiver view factor calculation method.
	- self.sim(): performs a combined simulation of both emissive and reflective losses and returns a global efficiency.
	- self.reset_opt(): empties the optics_manager bins to reuse the same object in a new ray trace.
	'''
	def __init__(self, apertureRadius, frustaRadii, frustaDepths, coneDepth, eps_wall, absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma, aperture_position, envelope_radius = None, envelope_depth = None, specular_receiver=False, dishDiameter_in = 20., sigma_in = 1.95e-3):
		'''
		Initialisation of the geometry and surface properties of the scene.

		The concentrator is a dish, the receiver is defined by a set of 2N parameters and housed in a cylindrical case.

		Arguments:
		- apertureRadius: radius of the aperture in (m).
		- frustaRadii: list of radii for the frusta elements in (m).
		- frustaDepths: list of depths for the frusta elements in (m).
		- coneDepth: depth of the cone element closing the geometry in (m).
		- eps_wall: minimum thickness of the insulation layer around the receiver in (m).
		- absReceiver: absorptivity of the walls of the cavity in the visible wavelength region. Can be a single value for a constant wall absorptivity or an array/list for a description of absorptivities element per element.
		- emsReceiver: emissivity of the walls of the cavity in the thermal emissions wavelength region. Can be a single value for a constant wall emissivity or an array/list for a description of emissivities element per element.
		- dishDiameter: diameter of the dish concentrator in (m).
		- dishFocus: focal distance of the dish concentrator in (m).
		- absDish: absorptivity of the dish mirrors.
		- sigma: surface slope error of the outer layer of the dish concentrator in (rad).
		- aperture_position: position of the receiver aperture plane on the z axis, qaxis of symmetry of the receiver and dish in (m).
		- envelope_radius: overrides the wall thickness based determination of the cylindrical envelope radius in (m).
		- envelope_depth: overrides the wall thickness based determination of the cylindrical envelope depth in (m).
		- specular_receiver: switches teh optics from lambertian to specular for the receiver internal walls.
		- dishDiameter_in: inner diameter of the dish in (m).
		- sigma_in: surface slope error of the inner layer of the dish concentrator in (rad).
		'''
		self.apertureRadius = apertureRadius
		self.frustaRadii = frustaRadii
		self.frustaDepths = frustaDepths # frusta depths
		self.coneDepth = coneDepth # cone depth
		self.dishDiameter = dishDiameter
		self.dishFocus = dishFocus
		self.absDish = absDish
		self.sigma = sigma
		if (type(absReceiver)!=N.ndarray) or (type(absReceiver)!=list):
			self.absReceiver = N.hstack((1,N.ones(len(frustaRadii)+1)*absReceiver))
		else:
			self.absReceiver = absReceiver
		assert(len(self.absReceiver)==(len(frustaRadii)+2))
		self.emsReceiver = emsReceiver
	
		self.Receiver = Assembly()
		Envelope = Assembly()
		self.Scene = Assembly()
		self.dishpos = -self.dishFocus
		self.aperture_position = aperture_position


		# Dish ----------------------------------------
		trd = translate(z=self.dishpos) # Dish frame transformation
		DISH = AssembledObject(surfs=[Surface(ParabolicDishGM(self.dishDiameter, self.dishFocus), RealReflectiveReceiver(self.absDish, self.sigma))], transform=trd)
		DISH2 = AssembledObject(surfs=[Surface(ParabolicDishGM(dishDiameter_in, self.dishFocus), RealReflectiveReceiver(self.absDish, sigma_in))], transform=translate(z=-dishFocus+0.0001))
		self.Scene.add_object(DISH)
		self.Scene.add_object(DISH2)

		# Envelope ----------------------------------------
		trr = translate(z=aperture_position) # Receiver frame transformation

		# Envelope maximum radius is the maximum between the aperture outer radius (to allow for pancakes) and the radii of frusta +the wall insulation thickness.
		if envelope_radius == None:		
			self.envelope_radius = N.amax(N.append(self.apertureRadius,N.array(self.frustaRadii)+eps_wall))
		else:
			self.envelope_radius = envelope_radius
		# Envelope maximum depth is the sum of all depths of the frusta + the wall insulation thickness. If the cone goes out, its depth has to be added too.
		if envelope_depth == None:
			max_depth = N.amax(N.add.accumulate(N.append(self.frustaDepths, self.coneDepth)))
			self.envelope_depth = max_depth+eps_wall
		else: 
			self.envelope_depth = envelope_depth

		if self.envelope_radius != self.apertureRadius:
			Envelope_front = AssembledObject(surfs=[Surface(RoundPlateGM(Ri=self.apertureRadius, Re=self.envelope_radius), LambertianReceiver(1))], transform=None)
			Envelope.add_object(Envelope_front)

		Envelope_back = AssembledObject(surfs=[Surface(RoundPlateGM(Re=self.envelope_radius), LambertianReceiver(1))], transform=translate(z=self.envelope_depth))
		Envelope_cylinder = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=2.*self.envelope_radius, height=self.envelope_depth), LambertianReceiver(1))], transform=translate(z=self.envelope_depth/2.))

		Envelope.add_object(Envelope_back)
		Envelope.add_object(Envelope_cylinder)

		# Receiver inner surface ----------------------------------------
		if specular_receiver==False:		
			# Cone
			if self.coneDepth < 0.: # self.coneDepth < 0 Inward cone
				CON = AssembledObject(surfs=[Surface(FiniteCone(r=self.frustaRadii[-1], h=-self.coneDepth), LambertianReceiver(self.absReceiver[-1]))], transform=translate(z=N.sum(frustaDepths)+self.coneDepth))
	
			elif self.coneDepth == 0.: # Round flat plates
				CON = AssembledObject(surfs=[Surface(RoundPlateGM(Re=self.frustaRadii[-1]), LambertianReceiver(self.absReceiver[-1]))], transform=translate(z=N.sum(frustaDepths)))

			else:	 # == cone depth > 0: Outgoing cone
				trc = N.dot(rotx(N.pi), translate(z=-N.sum(frustaDepths)-self.coneDepth)) # Cone frame transformation
				CON = AssembledObject(surfs=[Surface(FiniteCone(r=self.frustaRadii[-1], h=self.coneDepth), LambertianReceiver(self.absReceiver[-1]))], transform=trc)

			FRU = []
			# 1st frustum:
			if self.apertureRadius==self.frustaRadii[0]: # Cylinder
				frustum = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=self.apertureRadius*2, height=self.frustaDepths[0]), LambertianReceiver(self.absReceiver[1]))], transform=translate(z=self.frustaDepths[0]/2.))
			else:
				frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.apertureRadius, z2=self.frustaDepths[0], r2=self.frustaRadii[0]), LambertianReceiver(absorptivity=self.absReceiver[1]))], transform=None)
			FRU.append(frustum)
			# next frusta:
			for i in xrange(1,len(frustaRadii)):
				if self.frustaRadii[i-1]==self.frustaRadii[i]: # Cylinder
					frustum = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=self.frustaRadii[i-1]*2, height=self.frustaDepths[i]), LambertianReceiver(self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])+self.frustaDepths[i]/2.))
				elif self.frustaDepths[i] < 0.:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=-self.frustaDepths[i], r2=self.frustaRadii[i]), LambertianReceiver(absorptivity=self.absReceiver[1+i]))], transform=N.dot(translate(z=N.sum(self.frustaDepths[:i])),rotx(N.pi)))
				else:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=self.frustaDepths[i], r2=self.frustaRadii[i]), LambertianReceiver(absorptivity=self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])))
				FRU.append(frustum)

		if specular_receiver==True:
			# Cone
			if self.coneDepth < 0.: # self.coneDepth < 0 Inward cone
				CON = AssembledObject(surfs=[Surface(FiniteCone(r=self.frustaRadii[-1], h=-self.coneDepth), ReflectiveReceiver(self.absReceiver[-1]))], transform=translate(z=N.sum(frustaDepths)+self.coneDepth))

			elif self.coneDepth == 0.: # Round flat plates
				CON = AssembledObject(surfs=[Surface(RoundPlateGM(Re=self.frustaRadii[-1]), ReflectiveReceiver(self.absReceiver[-1]))], transform=translate(z=N.sum(frustaDepths)))

			else:	 # == cone depth > 0: Outgoing cone
				trc = N.dot(rotx(N.pi), translate(z=-N.sum(frustaDepths)-self.coneDepth)) # Cone frame transformation
				CON = AssembledObject(surfs=[Surface(FiniteCone(r=self.frustaRadii[-1], h=self.coneDepth), ReflectiveReceiver(self.absReceiver[-1]))], transform=trc)

			FRU = []
			# 1st frustum:
			if self.apertureRadius==self.frustaRadii[0]: # Cylinder
				frustum = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=self.apertureRadius*2, height=self.frustaDepths[0]), ReflectiveReceiver(self.absReceiver[1]))], transform=translate(z=self.frustaDepths[0]/2.))
			else:
				frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.apertureRadius, z2=self.frustaDepths[0], r2=self.frustaRadii[0]), ReflectiveReceiver(absorptivity=self.absReceiver[1]))], transform=None)
			FRU.append(frustum)
			# next frusta:
			for i in xrange(1,len(frustaRadii)):
				if self.frustaRadii[i-1]==self.frustaRadii[i]: # Cylinder
					frustum = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=self.frustaRadii[i-1]*2, height=self.frustaDepths[i]), ReflectiveReceiver(self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])+self.frustaDepths[i]/2.))
				elif self.frustaDepths[i] < 0.:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=-self.frustaDepths[i], r2=self.frustaRadii[i]), ReflectiveReceiver(absorptivity=self.absReceiver[1+i]))], transform=N.dot(translate(z=N.sum(self.frustaDepths[:i])),rotx(N.pi)))
				else:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=self.frustaDepths[i], r2=self.frustaRadii[i]), ReflectiveReceiver(absorptivity=self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])))
				FRU.append(frustum)

		# Assembly ----------------------------------------
		self.DISH = DISH
		self.DISH2 = DISH2
		self.CON = CON
		self.FRU = FRU
		self.Receiver.add_object(CON)	
		for i in xrange(len(FRU)):
			self.Receiver.add_object(self.FRU[i])
		self.Scene.add_assembly(Envelope, transform=trr)
		self.Scene.add_assembly(self.Receiver, transform=trr)

	def ray_sim(self, nrays, G=1000., CSR=None):
		'''
		Raytracing-only method.

		Arguments:
		- nrays: number of rays in the bundle
		- G: DNI in (W/m2)
		- CSR: Circumsolar ratio for the Buie sunshape. If it stays as None, a pillbox sunshape model is usend instead.

		Returns:
		- An array including: the 2N parameters in (m), flux input and losses breakdown in (W).
		The hits and absorbed energy are stored in the following variables:
			- self.Frusta_hits, self.Frusta_abs: location and value of the absorbed energy at this location on the frusta.
			- self.Cone_hits, self.Cone_abs: location and value of the absorbed energy at this location on the cone.
		'''
		self.nrays = nrays
		
		# Source declaration
		sourceCenter = N.array([[0,0,2.*self.dishFocus]]).T # Source center position
		sourceDirection = N.array([0,0,-1.]) # Source normal direction
		sourceRadius = 0.6*self.dishDiameter # m, Source radius
		sourceAngle = 4.65e-3 # radians, sun rays angular range
		if CSR == None: # choice between Buie or pillbox sunshapes.
			SOURCE = solar_disk_bundle(nrays, sourceCenter, sourceDirection, sourceRadius, sourceAngle, G)
		else:
			self.CSR = CSR
			SOURCE = buie_sunshape(nrays, sourceCenter, sourceDirection, sourceRadius, CSR, G)

		# RAYTRACE!
		self.engine = TracerEngine(self.Scene)
		itmax = 30 # stop iteration after this many ray bundles were generated (i.e. 
					# after the original rays intersected some surface this many times).
		self.minener = 1e-10 # minimum energy threshold
		self.engine.ray_tracer(SOURCE, itmax, self.minener, tree=True)
		
		# RESULTS:
		# Absorbed flux
		Dish_abs, Dish_hits = self.DISH.get_surfaces()[0].get_optics_manager().get_all_hits()
		Dish2_abs, Dish2_hits = self.DISH2.get_surfaces()[0].get_optics_manager().get_all_hits()
		Cone_abs, Cone_hits = self.CON.get_surfaces()[0].get_optics_manager().get_all_hits()
		Frusta_abs = []
		Frusta_abs_total = 0
		Frusta_hits = []
		for i in xrange(len(self.FRU)):
			Frustum_abs, Frustum_hits = self.FRU[i].get_surfaces()[0].get_optics_manager().get_all_hits()
			Frusta_abs_total += N.sum(Frustum_abs)
			Frusta_abs.append(Frustum_abs)
			Frusta_hits.append(Frustum_hits)
		
		self.ReceiverA = N.sum(Cone_abs)+Frusta_abs_total
		self.Frusta_abs = Frusta_abs
		self.Frusta_hits = Frusta_hits
		self.Cone_abs = Cone_abs
		self.Cone_hits = Cone_hits
		
		# Blockage and spillage losses:
		
		# Flux input (just considering useful rays from the source).
		self.flux_input = N.sum(self.engine.tree._bunds[0].get_energy()[self.engine.tree._bunds[1].get_parents()])

		# Receiver reflectivity losses:
		# Is considered receiver reflectivity losses all that is lost from the raytrace and not identified as one of the previous losses.
		rec_ref = 0
		spill = 0	
		block = 0

		for i in xrange(1,len(self.engine.tree._bunds)-1):
			bund_i_rays = self.engine.tree._bunds[i].get_num_rays()
			bund_i_parents = self.engine.tree._bunds[i].get_parents()
			bund_i_ener = self.engine.tree._bunds[i].get_energy()
			bund_im1_ener = self.engine.tree._bunds[i-1].get_energy()
			bund_ip1_parents = self.engine.tree._bunds[i+1].get_parents()

			hits = N.zeros(bund_i_rays)
			# hits = 1 means that actual rays hit something and will exist in the next bundle.
			hits[bund_ip1_parents] = 1
			# miss = 1 means that the ray is not hitting anything and will not exist in the next bundle.
			miss = hits != 1
			# disappear means that the ray has hit something and was completely lost afterwards, as a consequence we know it didn't hit an active surface.
			disappeared = bund_i_parents[bund_i_ener < self.minener]

			if i==1:
				block += N.sum(bund_im1_ener[disappeared])
				spill += N.sum(bund_i_ener[miss])
			elif i==2:
				spill += N.sum(bund_im1_ener[disappeared])
				rec_ref += N.sum(bund_i_ener[miss])
			else:
				rec_ref += N.sum(bund_i_ener[miss])
		
		# Dish absorptivity losses: small imprecision here due to blocked rays absorbed by the dish in bundle 2.
		self.dish_abs_losses = N.sum(Dish_abs)+N.sum(Dish2_abs) # all that hits directly from the source and loses energy is due to dish absorptivity losses.

		self.block_losses = block
		self.rec_ref_losses = rec_ref
		self.spill_losses = spill

		# Optical losses total:
		self.opt_losses_total = self.flux_input-self.ReceiverA

		return [self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, self.flux_input, self.ReceiverA, self.opt_losses_total ,self.block_losses, self.spill_losses, self.dish_abs_losses, self.rec_ref_losses]

	def temperature_guess(self, T_in, p_in, T_out, tube_diameters_in, tube_diameters_out, tube_conductivity, emissions_guess, passive = None, coating_thickness = 45e-6, coating_conductivity = 1.2):
		'''
		Makes a first guess on temperature profiles approximating the enthalpy gain of the water/steam mixture to be equal to flux input on the tubes external walls. The tube walls are coated with a selective coating. Default arguments are for Pyromark2500(R).

		Arguments:
		- T_in: Inlet temperature of the water  in (K).
		- p_in: Inlet pressure of the water in (pa).
		- T_out: Outlet temeprature of the water in (K).
		- tube_diameters_in: inner diameter of the tubes in (m).
		- tube_diameters_out: outer diameter of the tubesin (m).
		- tube_conductivity: thermal conductivity of teh tubes in (W/mK).
		- emissions_guess: emissive losses guess to account for thermal emissions in (W).
		- passive: array of the indices of the adiabatic surfaces in the cavity.
		- coating_thickness: thickness of the coating layer on the tubes in (m).
		- coating_conductivity: coating thermal conductivity in (W/mK).

		Returns:
		- strings 'good_geom' or 'bad_geom' depending on the mass flow guess to meet the input/output arguments and amount of actually going in the receiver. This is a quick hack to prevent issues with receivers forcing in the required inpu/output by lowering the mass flow too much/putting it negative, thus impacting the enthalpy guess... and basically screwing-up the convergence process for non-performing geometries.
		'''
		# Get starting and ending enthalpies via Freesteam
		h_in = steam_pT(p_in,T_in).h
		p_out = 0.9*p_in # arbitrary 10% of pressure drop based on first 4 parameters cavity cases observation

		# Detect active surfaces
		active = N.ones(len(self.areas), dtype = N.bool)
		active[0] = 0 # aperture
		if passive != None:
			active[passive] = 0

		# Evaluate pressure drop per tube/sections
		tube_lengths = N.array(self.areas[active]/tube_diameters_out)
		ps = [p_in]
		ps.append(p_in+(p_out-p_in)/N.sum(tube_lengths)*N.add.accumulate(tube_lengths))
		ps = N.hstack(ps)

		# Build a bin matching enthalpy array
		hs = N.zeros(len(self.areas[active])+1)
		self.rad_bins = N.zeros(len(self.areas)-1)
		hs[0] = h_in
		index = 0 # to track which element of the profile we are dealing with in absolute terms.
		active_abs = N.zeros(len(self.areas)) # To track the radiative energy absorbed by the active zone in total.
		# Bin Frusta:
		for i in xrange(len(self.frustaDepths)):
			z1 = N.sum(self.frustaDepths[:i+1])-self.frustaDepths[i]
			z2 = N.sum(self.frustaDepths[:i+1])
			# Adjust the hits detector position to take into account the placement of the receicer in comparision with the focal plane.
			z1 += self.aperture_position
			z2 += self.aperture_position
			if self.frustaDepths[i]<0.:
					z1,z2 = z2,z1

			for j in xrange(self.bins_frustum[i]):
				index+=1
				section = N.logical_and(self.Frusta_hits[i][2]>(z1+(z2-z1)*j/self.bins_frustum[i]), self.Frusta_hits[i][2]<=(z1+(z2-z1)*(j+1)/self.bins_frustum[i]))
				self.rad_bins[index-1] = N.sum(self.Frusta_abs[i][section])
				if active[index]:
					active_abs[index] = N.sum(self.Frusta_abs[i][section])
		# Bin the cone:
		for i in xrange(self.bins_cone):
			index+=1
			r1 = self.frustaRadii[-1]-i*self.frustaRadii[-1]/self.bins_cone
			r2 = self.frustaRadii[-1]-(i+1)*self.frustaRadii[-1]/self.bins_cone
			cone_hits_radii = N.sqrt(self.Cone_hits[0]**2+self.Cone_hits[1]**2)
			section = N.logical_and(cone_hits_radii<r1, cone_hits_radii>=r2)
			self.rad_bins[index-1] = N.sum(self.Cone_abs[section])
			if active[index]:
				active_abs[index] = N.sum(self.Cone_abs[section])

		# Evaluate the mass flow using active zones only and subtracting the emissive losses:
		h_out = steam_pT(p_out,T_out).h

		emissions_active = N.zeros(len(active_abs))
		emissions_active[active] = emissions_guess[active]
		self.m = N.sum(active_abs[active]+emissions_guess[active])/(h_out-h_in)
		#print 'm =', self.m
		#print 'rad_bins', self.rad_bins, 'total abs',  N.sum(self.rad_bins)

		if self.m < 0.01:
			return 'bad_geom'

		# Evlauate the enthalpies
		for i in xrange(1,len(hs)):
			hs[i] = hs[i-1]+(active_abs[i]+emissions_active[i])/self.m

		# Get temperatures from enthalpies via Freesteam
		self.T_guess_fluid = N.zeros(len(self.areas))
		self.T_guess_fluid[0] = T_in
		for i in xrange(1,len(hs)):
			self.T_guess_fluid[i] = steam_ph(ps[i],hs[i]).T

		#print 'emissions_guess', emissions_guess
		T_guess_wall = self.T_guess_fluid
		T_guess_wall[active]+=(active_abs[active]+emissions_guess[active])/(N.pi*tube_lengths)*(N.log(tube_diameters_out/tube_diameters_in)/tube_conductivity+N.log((tube_diameters_out+2.*coating_thickness)/tube_diameters_out)/coating_conductivity)
		T_guess_wall[~active] = T_out
		T_guess_wall[0] = T_in
		#print 'T wall:', T_guess_wall
		self.T_guess = ((T_guess_wall[1:]**4.+T_guess_wall[:-1]**4.)/2.)**(1./4.)
		#print 'T final', self.T_guess
		assert (self.T_guess==float('inf')).any()==False, str(self.T_guess)+str(T_guess_wall)+str(emissions_guess)+str(self.m)+str([self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth])

		self.rad_passive = N.zeros(N.shape(self.rad_bins))
		if passive != None:
			self.rad_passive[passive] = self.rad_bins[passive]

		return 'good_geom'

	def VF_sim(self, bins_frustum, bins_cone, precision=0.005):
		'''
		Calculates view factors and areas for a Two_N_parameters_cavity instance.
		
		Arguments:
		- bins_frustum: array/list of the number of elements per frusta sections used in the radiative power binning process. Every bin on a single frustum element has the same depth on the z axis (symmetry axis of the receiver geometry).
		- bins_cone: number of elements per frusta sections used in the radiative power binning process. Every bin on the cone section has the same depth on the z axis (symmetry axis of the receiver geometry).
		- precision: precision required for the view factor calculations.

		Returns:
		- self.VF: view factor matrix.
		- self.areas: array containing the areas of each of the bins in (m2).
		'''
		# Emissivity array shaping. Aperture value has to be 1 for an open receiver. Takes place before the VF calculation to avoid spending useless time on this if the array is not shaped properly.
		if (type(self.emsReceiver)!=N.ndarray) and (type(self.emsReceiver)!=list):
			self.emsReceiver = N.hstack(N.append(1,N.ones(N.sum(N.array(bins_frustum))+bins_cone)*self.emsReceiver))
		assert(len(self.emsReceiver)==N.sum(N.array(bins_frustum))+bins_cone+1)

		# Calculate geometry view factor
		vfcase = Two_N_parameters_cavity_RTVF(self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, el_FRUs=bins_frustum, el_CON=bins_cone, precision=precision)

		self.bins_frustum = bins_frustum
		self.bins_cone = bins_cone
		self.VF = vfcase.VF_esperance
		self.areas = vfcase.areas
		
		return self.VF, self.areas

	def emi_sim(self, Tamb, Trec, VF, areas, inc_radiation=None):
		'''
		Method to calculate the emissive losses of a Two_N_parameters_cavity instance. Calls a radiosity method using the view factors matrix calculated previously

		Arguments:
		- Tamb: temperature of the environment in (K), used for the aperture in the radiosity method.
		- Trec: array/list of temperatures of the surfaces of the receiver in (K). If only one value is declared, set all receiver surfaces to this single temperature value.
		- VF: View factors matrix.
		- areas: Areas of the elemenst of the view factor matrix (m2).
		- inc_radiation: array/list of incoming radiative power on the elements (W). If not None overrides the radiosity problem temperature boundary condition only where the inc_radiation value is not equal to 0.

		Returns: 
		- self.emissive_losses: overall emissive losses through the aperture
		The method also stores the temperatures and net emitted radiative flux and net radiative powers in the thermal wavelength region (semi-gray body assumption):
			- self.q: Net thermal radiative flux (W/m2) per element.
			- self.Q: Net thermal radiative power (W) per element.
			- self.T: Temeprature of each element (K).
		'''
		self.tempAmbiant = Tamb
		self.tempReceiver = Trec

		if inc_radiation != None:
			inc_radiation = N.hstack((0.,inc_radiation)) # To take into account the aperture in the radiosity system.
		# Solve radiosity problem
		AA,bb,J,Eb,T,q,Q = radiosity_RTVF(VF, areas, self.emsReceiver, Tamb, Trec, inc_radiation)
		self.q = q
		self.Q = Q
		self.T = T		

		self.emissive_losses = -self.Q[0]

		return self.emissive_losses

	def reset_opt(self):
		S = self.engine._asm.get_surfaces()
		for i in xrange(len(S)):
			S[i].get_optics_manager().reset()

	def sim(self, Tamb, Trec, nrays, CSR=None):
		'''
		Method to simulate the radiative efficiency of a Two_N_parameters_cavity receiver with fixed temperatures boundary condition, including ray tracing of visible light and thermal emeissions.

		Arguments:
		- Tamb: temperature of the environment in (K), used for the aperture in the radiosity method.
		- Trec: array/list of temperatures of the surfaces of the receiver in (K). If only one value is declared, set all receiver surfaces to this single temperature value.
		- nrays: number of rays in the bundle
		- CSR: Circumsolar ratio for the Buie sunshape. If it stays as None, a pillbox sunshape model is usend instead.

		Returns:
		- An array including: the 2N parameters in (m), flux input and losses breakdown in (W).
		The hits and absorbed energy are stored in the following variables:
			- self.Frusta_hits, self.Frusta_abs: location and value of the absorbed energy at this location on the frusta.
			- self.Cone_hits, self.Cone_abs: location and value of the absorbed energy at this location on the cone.
		The method also stores the temperatures and net emitted radiative flux and net radiative powers in the thermal wavelength region (semi-gray body assumption):
			- self.q: Net thermal radiative flux (W/m2) per element.
			- self.Q: Net thermal radiative power (W) per element.
			- self.T: Temeprature of each element (K).
		'''

		self.ray_sim(nrays, CSR=CSR)
		self.emi_sim(Tamb, Trec, VF=self.VF, areas=self.areas, inc_radiation=None)

		System_efficiency = (self.flux_input-self.opt_losses_total-self.emissive_losses)/self.flux_input

		return [self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, self.flux_input, self.ReceiverA, self.emissive_losses, self.opt_losses_total ,self.block_losses, self.spill_losses, self.dish_abs_losses, self.rec_ref_losses]

	def sim_with_T_guess(self, Tamb, Trec_in, p_in, Trec_out, tube_diameters_in, tube_diameters_out, tube_conductivity, nrays, G=1000., CSR=None, passive = None):
		'''
		Method to simulate the radiative efficiency of a Two_N_parameters_cavity receiver with a realistic evaluation of the temepratures of the walls using fluid properties and the heat exchange model from the temperature_guess() method.

		Arguments:
		- T_in: Inlet temperature of the water  in (K).
		- p_in: Inlet pressure of the water in (pa).
		- T_out: Outlet temeprature of the water in (K).
		- tube_diameters_in: inner diameter of the tubes in (m).
		- tube_diameters_out: outer diameter of the tubesin (m).
		- tube_conductivity: thermal conductivity of teh tubes in (W/mK).
		- nrays: number of rays in the bundle
		- G: DNI in (W/m2)
		- CSR: Circumsolar ratio for the Buie sunshape. If it stays as None, a pillbox sunshape model is usend instead.
		- passive: array of the indices of the adiabatic surfaces in the cavity.

		Returns:
		- An array including: the 2N parameters in (m), flux input and losses breakdown in (W).
		The hits and absorbed energy are stored in the following variables:
			- self.Frusta_hits, self.Frusta_abs: location and value of the absorbed energy at this location on the frusta.
			- self.Cone_hits, self.Cone_abs: location and value of the absorbed energy at this location on the cone.
		The method also stores the temperatures and net emitted radiative flux and net radiative powers in the thermal wavelength region (semi-gray body assumption):
			- self.q: Net thermal radiative flux (W/m2) per element.
			- self.Q: Net thermal radiative power (W) per element.
			- self.T: Temeprature of each element (K).
		'''
		# Optical ray trace for visible light
		self.ray_sim(nrays, CSR=CSR, G=G)

		# Iterate to find wall temperatures and respect the energy balance:
		if type(self.emsReceiver)==float:
			self.emsReceiver = N.hstack((1,N.ones(len(self.areas)-1)*self.emsReceiver))
		emissions = N.ones(len(self.areas))
		convergence = N.ones(len(emissions))
		while (convergence>0.001).any():
			result_T_guess = self.temperature_guess(Trec_in, p_in, Trec_out, tube_diameters_in, tube_diameters_out, tube_conductivity, emissions, passive)
			if result_T_guess == 'bad_geom': # discard 'bad_geom' geometries.
				[self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, self.flux_input, self.ReceiverA, self.emissive_losses, self.opt_losses_total ,self.block_losses, self.spill_losses, self.dish_abs_losses, self.rec_ref_losses] = [self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, self.flux_input, 0., 0., 0. ,0., 0., 0., 0.]
				break
			self.emi_sim(Tamb, self.T_guess, VF=self.VF, areas=self.areas, inc_radiation=self.rad_passive)
			self.T_guess = self.T
			#print 'Final T guess:', self.T_guess
			convergence = N.abs((self.Q-emissions)/self.Q)
			emissions = (self.Q+emissions)/2.
			#print '     '
		self.reset_opt()

		return [self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, self.flux_input, self.ReceiverA, self.emissive_losses, self.opt_losses_total ,self.block_losses, self.spill_losses, self.dish_abs_losses, self.rec_ref_losses]
