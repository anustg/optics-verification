'''
__________________________________________________________________________________________________
*********** 2 N Parameters cavity receiver  ************
__________________________________________________________________________________________________
'''
import numpy as N

from tracer.surface import *
from tracer.quadric import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *

from freesteam import *

from emissive_losses.emissive_losses import *
from emissive_losses.view_factors_3D import *

class TwoNparamcav(Assembly):
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
	def __init__(self, apertureRadius, frustaRadii, frustaDepths, coneDepth, eps_wall, absReceiver, emsReceiver, aperture_position, envelope_radius = None, envelope_depth = None, specular_receiver=False):
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

		- aperture_position: position of the receiver aperture plane on the z axis, qaxis of symmetry of the receiver and dish in (m).
		- envelope_radius: overrides the wall thickness based determination of the cylindrical envelope radius in (m).
		- envelope_depth: overrides the wall thickness based determination of the cylindrical envelope depth in (m).
		- specular_receiver: switches teh optics from lambertian to specular for the receiver internal walls.

		'''
		Assembly.__init__(self, objects=None, subassemblies=None, location=None, rotation=None)

		self.apertureRadius = apertureRadius
		self.frustaRadii = frustaRadii
		self.frustaDepths = frustaDepths # frusta depths
		self.coneDepth = coneDepth # cone depth

		if (type(absReceiver)==float) or (type(absReceiver)==int):
			self.absReceiver = N.hstack((1.,N.ones(len(frustaRadii)+1)*absReceiver))
		else:
			self.absReceiver = absReceiver
		assert(len(self.absReceiver)==(len(frustaRadii)+2))
		self.emsReceiver = emsReceiver
	
		Envelope = Assembly()
		self.Active_zone = Assembly()

		self.aperture_position = aperture_position

		# Envelope ----------------------------------------
		trr = translate(z=aperture_position) # Receiver frame transformation

		# Envelope maximum radius is the maximum between the aperture outer radius (to allow for pancakes) and the radii of frusta +the wall insulation thickness.
		if envelope_radius == None:		
			self.envelope_radius = N.amax(N.hstack([self.apertureRadius,N.array(self.frustaRadii)+eps_wall]))
		else:
			self.envelope_radius = envelope_radius
		# Envelope maximum depth is the sum of all depths of the frusta + the wall insulation thickness. If the cone goes out, its depth has to be added too.
		if envelope_depth == None:
			max_depth = N.amax(N.add.accumulate(N.hstack([self.frustaDepths, self.coneDepth])))
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
			elif self.frustaDepths[0] == 0.:
				frustum = AssembledObject(surfs=[Surface(RoundPlateGM(Re=self.apertureRadius, Ri=self.frustaRadii[0]), LambertianReceiver(self.absReceiver[1]))], transform=translate(z=N.sum(self.frustaDepths[:i])))
			else:
				frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.apertureRadius, z2=self.frustaDepths[0], r2=self.frustaRadii[0]), LambertianReceiver(absorptivity=self.absReceiver[1]))], transform=None)
			FRU.append(frustum)
			# next frusta:
			for i in xrange(1,len(frustaRadii)):
				if self.frustaRadii[i-1]==self.frustaRadii[i]: # Cylinder
					frustum = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=self.frustaRadii[i-1]*2, height=self.frustaDepths[i]), LambertianReceiver(self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])+self.frustaDepths[i]/2.))
				elif self.frustaDepths[i] < 0.:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=-self.frustaDepths[i], r2=self.frustaRadii[i]), LambertianReceiver(absorptivity=self.absReceiver[1+i]))], transform=N.dot(translate(z=N.sum(self.frustaDepths[:i])),rotx(N.pi)))
				elif self.frustaDepths[i] > 0.:
					frustum = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0., r1=self.frustaRadii[i-1], z2=self.frustaDepths[i], r2=self.frustaRadii[i]), LambertianReceiver(absorptivity=self.absReceiver[1+i]))], transform=translate(z=N.sum(self.frustaDepths[:i])))
				else:
					frustum = AssembledObject(surfs=[Surface(RoundPlateGM(Re=self.frustaRadii[i-1], Ri=self.frustaRadii[i]), ReflectiveReceiver(self.absReceiver[-1]))], transform=translate(z=N.sum(self.frustaDepths[:i])))

				FRU.append(frustum)
				
		# Receiver Assembly ----------------------------------------

		self.CON = CON
		self.FRU = FRU
		for i in xrange(len(FRU)):
			self.Active_zone.add_object(FRU[i])
		self.Active_zone.add_object(CON)	
		self.add_assembly(Envelope, transform=trr)
		self.add_assembly(self.Active_zone, transform=trr)

	def VF_sim(self, bins_frusta, bins_cone, num_rays, precision):
		'''
		Calculates view factors and areas for a Two_N_parameters_cavity instance.
		
		Arguments:
		- bins_frusta: array/list of the number of elements per frusta sections used in the radiative power binning process. Every bin on a single frustum element has the same depth on the z axis (symmetry axis of the receiver geometry).
		- bins_cone: number of elements per frusta sections used in the radiative power binning process. Every bin on the cone section has the same depth on the z axis (symmetry axis of the receiver geometry).
		- precision: precision required for the view factor calculations.

		Returns:
		- self.VF: view factor matrix.
		- self.areas: array containing the areas of each of the bins in (m2).
		'''
		# Emissivity array shaping. Aperture value has to be 1 for an open receiver. Takes place before the VF calculation to avoid spending useless time on this if the array is not shaped properly.
		if (type(self.emsReceiver)!=N.ndarray) and (type(self.emsReceiver)!=list):
			self.emsReceiver = N.hstack(N.append(1.,N.ones(N.sum(N.array(bins_frusta))+bins_cone)*self.emsReceiver))
		assert(len(self.emsReceiver)==N.sum(N.array(bins_frusta))+bins_cone+1)

		# Calculate geometry view factor
		vfcase = Two_N_parameters_cavity_RTVF(self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth, el_FRUs=bins_frusta, el_CON=bins_cone, num_rays=num_rays, precision=precision)

		self.bins_frusta = bins_frusta
		self.bins_cone = bins_cone
		self.VF = vfcase.VF_esperance
		self.areas = vfcase.areas

		return self.VF, self.areas, self.bins_frusta, self.bins_cone

	def bin_hits(self):

		int_walls = self.Active_zone.get_surfaces()
		self.bin_abs = N.zeros(len(self.areas)-1)
		index = 0 # to track which element of the profile we are dealing with in absolute terms.
		active_abs = N.zeros(len(self.areas)) # To track the radiative energy absorbed by the active zone in total.

		receiver_abs = []
		receiver_hits = []

		# Bin Frusta:
		for i in xrange(len(self.frustaDepths)):
			z1 = N.sum(self.frustaDepths[:i+1])-self.frustaDepths[i]
			z2 = N.sum(self.frustaDepths[:i+1])
			if i == 0:
				r1 = self.apertureRadius
			else:
				r1 = self.frustaRadii[i-1]
			r2 = self.frustaRadii[i]
			# Adjust the hits detector position to take into account the placement of the receiver in comparision with the focal plane.
			z1 += self.aperture_position
			z2 += self.aperture_position
			if self.frustaDepths[i]<0.:
				z1,z2 = z2,z1
			if r1>r2:
				r1,r2 = r2,r1

			gethits = int_walls[i].get_optics_manager().get_all_hits()
			abs = gethits[0]
			hits = gethits[1]
			receiver_abs.append(abs)
			receiver_hits.append(hits)
			
			for j in xrange(self.bins_frusta[i]):
				index+=1
				z1bin = z1+(z2-z1)*j/self.bins_frusta[i]
				z2bin = z1+(z2-z1)*(j+1)/self.bins_frusta[i]

				test_depth = N.logical_and(hits[2]>z1bin, hits[2]<=z2bin)

				self.bin_abs[index-1] = N.sum(abs[test_depth])

		# Bin the cone:
		gethits = int_walls[-1].get_optics_manager().get_all_hits()
		abs = gethits[0]
		hits = gethits[1]
		receiver_abs.append(abs)
		receiver_hits.append(hits)
		for i in xrange(self.bins_cone):
			index+=1
			r1 = self.frustaRadii[-1]-i*self.frustaRadii[-1]/self.bins_cone
			r2 = self.frustaRadii[-1]-(i+1)*self.frustaRadii[-1]/self.bins_cone
			cone_hits_radii = N.sqrt(hits[0]**2+hits[1]**2)
			section = N.logical_and(cone_hits_radii<r1, cone_hits_radii>=r2)
			self.bin_abs[index-1] = N.sum(abs[section])

		receiver_abs = N.hstack(receiver_abs)
		receiver_hits = N.concatenate(receiver_hits, axis=1)		

		return self.bin_abs, receiver_abs, receiver_hits


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
		active[0] = 0. # aperture
		if passive != None:
			active[passive] = 0.
		active_abs = N.zeros(len(self.areas))
		active_abs[1:] = self.bin_abs*active[1:]

		# Evaluate pressure drop per tube/sections
		tube_lengths = N.array(self.areas[active]/tube_diameters_out)
		self.tube_lengths = tube_lengths
		ps = [p_in]
		ps.append(p_in+(p_out-p_in)/N.sum(tube_lengths)*N.add.accumulate(tube_lengths))
		ps = N.hstack(ps)

		# Build a bin matching enthalpy array
		hs = N.zeros(len(self.areas[active])+1)
		hs[0] = h_in

		# Evaluate the mass flow using active zones only and subtracting the emissive losses:
		h_out = steam_pT(p_out,T_out).h

		emissions_active = N.zeros(len(self.areas))
		emissions_active[active] = emissions_guess[active]
		self.m = N.sum(active_abs[active]+emissions_guess[active])/(h_out-h_in)

		# Evlauate the enthalpies
		for i in xrange(1,len(hs)):
			hs[i] = hs[i-1]+(active_abs[i]+emissions_active[i])/self.m

		self.h = hs

		#FIXME: need a more reliable convergence insurance
		if self.m < 0.01:
			return 'bad_geom'

		# Get temperatures from enthalpies via Freesteam
		self.T_guess_fluid = N.zeros(len(self.areas))
		self.T_guess_fluid[0] = T_in
		for i in xrange(1,len(hs)):
			self.T_guess_fluid[i] = steam_ph(ps[i],hs[i]).T
		#print 'T fluid:', self.T_guess_fluid
		T_guess_wall = self.T_guess_fluid
		T_guess_wall[active]+=(active_abs[active]+emissions_guess[active])/(N.pi*tube_lengths)*(N.log(tube_diameters_out/tube_diameters_in)/tube_conductivity+N.log((tube_diameters_out+2.*coating_thickness)/tube_diameters_out)/coating_conductivity)
		T_guess_wall[~active] = T_out
		T_guess_wall[0] = T_in
		#print 'T wall:', T_guess_wall
		self.T_guess = ((T_guess_wall[1:]**4.+T_guess_wall[:-1]**4.)/2.)**(1./4.)
		#print 'T final', self.T_guess
		assert (self.T_guess==float('inf')).any()==False, str(self.T_guess)+str(T_guess_wall)+str(emissions_guess)+str(self.m)+str([self.apertureRadius, self.frustaRadii, self.frustaDepths, self.coneDepth])

		self.rad_passive = None
		if passive != None:
			self.rad_passive = N.zeros(N.shape(self.bin_abs))
			self.rad_passive.fill(N.nan)
			self.rad_passive[passive] = self.bin_abs[passive]

		return 'good_geom'

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
		if inc_radiation != None:
			inc_radiation = N.hstack((N.nan,inc_radiation)) # To take into account the aperture in the radiosity system.
		# Solve radiosity problem
		T = N.hstack((Tamb, Trec))
		AA,bb,J,Eb,T,q,Q = radiosity_RTVF(VF, areas, self.emsReceiver, T, inc_radiation)
		self.q = q
		self.Q = Q
		self.T = T		

		self.emissive_losses = -self.Q[0]

		return self.emissive_losses

	def energy_balance(self, Tamb, Trec_in, p_in, Trec_out, tube_diameters_in, tube_diameters_out, tube_conductivity, passive = None):
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

		- passive: array of the indices of the adiabatic surfaces in the cavity.

		Returns:
		The temperature of the elements, an array of zeros if the candidates are net energy destructors.
		'''

		# Iterate to find wall temperatures and respect the energy balance:
		if type(self.emsReceiver)==float:
			self.emsReceiver = N.hstack((1.,N.ones(len(self.areas)-1)*self.emsReceiver))

		emissions = N.ones(len(self.areas))
		convergence = N.ones(len(emissions))

		while (convergence>0.001).any():
			result_T_guess = self.temperature_guess(Trec_in, p_in, Trec_out, tube_diameters_in, tube_diameters_out, tube_conductivity, emissions, passive)

			if result_T_guess == 'bad_geom': # discard 'bad_geom' geometries.

				self.T_guess = N.ones(len(self.areas))*Trec_in
				break

			self.emi_sim(Tamb, self.T_guess, VF=self.VF, areas=self.areas, inc_radiation=self.rad_passive)
			self.T_guess = self.T
			#print self.T
			#print 'Final T guess:', self.T_guess
			convergence = N.abs((self.Q-emissions)/self.Q)
			emissions = (self.Q+emissions)/2.
		return result_T_guess
