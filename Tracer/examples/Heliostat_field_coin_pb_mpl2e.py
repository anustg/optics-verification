'''
A tower/heliostat field example based on the Tower_gui.py example.
It has been cleaned from GUI and mayaVi stuff and given the CoIn rendering import and capacity.
'''
import traits.api as t_api
import traitsui.api as tui
from tracer.CoIn_rendering.rendering import *

from tracer.ray_bundle import RayBundle, concatenate_rays
from tracer.sources import solar_disk_bundle, buie_sunshape
from tracer.assembly import Assembly
from tracer.spatial_geometry import rotx, roty, rotz, rotation_to_z
from tracer.tracer_engine import TracerEngine

from tracer.models.one_sided_mirror import one_sided_receiver
from tracer.models.heliostat_field import HeliostatField, radial_stagger, solar_vector

import numpy as N
from scipy.constants import degree
import math
import matplotlib.pyplot as plt


class TowerScene():

    # Location of the sun:
    sun_az = 80
    sun_elev = 45

    sun_vec = solar_vector(sun_az*degree, sun_elev*degree)
    hstat_normals = N.zeros((218,3))
    
    # import custom coordinate file
    pos = N.loadtxt("sandia_hstat_coordinates.csv", delimiter=',')

    # Field-based calculations for sunshape parameters
    #===================================================
    t_pos = pos.T
    xc_min = t_pos[0][N.argmin(t_pos[0])]
    xc_max = t_pos[0][N.argmax(t_pos[0])]
    yc_min = t_pos[1][N.argmin(t_pos[1])]
    yc_max = t_pos[1][N.argmax(t_pos[1])]

    x_dist = xc_max - xc_min
    y_dist = yc_max - yc_min

    xc_cent = (xc_min + xc_max) / 2
    yc_cent = (yc_min + yc_max) / 2
    field_centre = N.r_[xc_cent, yc_cent, 0]
    print('field centre', field_centre)
    #===================================================
    
    def __init__(self):
        self.gen_plant() 
   
    def gen_rays(self, num_rays):

        #========================
        individual_source = False
        #========================

        if individual_source:
            # Pillbox source on a per-heliostat basis
            radius = 1.20 * math.sqrt(2 * 3.405**2)

            direction = N.array(-self.sun_vec)

            ray_list = []
            num_surfs = self.pos.shape[0]
            for i in xrange(num_surfs):
                centre = N.c_[50 * self.sun_vec + self.pos[i]]
                rayb = solar_disk_bundle(num_rays/num_surfs, centre, direction, radius, 1e-3, 1000)
                ray_list.append(rayb)

            rays = concatenate_rays(ray_list)
            del ray_list
            
        else:
            # Giant source covering entire field area
            radius = 1.10 * math.sqrt((self.x_dist/2)**2 + (self.y_dist/2)**2)
            print('Radius', radius)
            global sun_area
            sun_area = N.pi * radius**2
            print('Size of source', sun_area)

            centre = N.c_[300*self.sun_vec + self.field_centre]
            direction = N.array(-self.sun_vec)

            rays = solar_disk_bundle(num_rays, centre, direction, radius, 4.65e-3, 1000)
                
        return rays
    
    def gen_plant(self):

        # set heliostat field characteristics: 6.09m*6.09m, abs = 0, aim_h = 61
        self.field = HeliostatField(self.pos, 6.09e-0, 6.09e-0, 0, 61)
        
        self.rec, recobj = one_sided_receiver(11, 11)
        rec_trans = rotx(N.pi/-2)
        rec_trans[2,3] = 61

        #=================
        ground_rec = False
        #=================

        if ground_rec:
            # Evaluating missed rays in the field along with receiver
            radius = 1.10 * math.sqrt((self.x_dist/2)**2 + (self.y_dist/2)**2)
            self.ground_rec, ground_recobj = one_sided_receiver(3*radius, 3*radius)

            ground_rec_trans = rotz(0)
            ground_rec_trans[0,3] = self.field_centre[0]
            ground_rec_trans[1,3] = self.field_centre[1]
            
            recobj.set_transform(rec_trans)
            ground_recobj.set_transform(ground_rec_trans)
            self.plant = Assembly(objects=[recobj, ground_recobj], subassemblies=[self.field])
        else:
            # Evaluating just the receiver
            recobj.set_transform(rec_trans)
        
            self.plant = Assembly(objects=[recobj], subassemblies=[self.field])      
    
    def aim_field(self):
        hstat_az, hstat_elev = \
            self.field.aim_to_sun(self.sun_az*degree, self.sun_elev*degree)
        print('Heliostat azimuths', (N.around(hstat_az, decimals=3)))
        print('Heliostat elevations', (N.around(hstat_elev, decimals=3)))

        return hstat_az, hstat_elev

    def calculate_area(self, hstat_az, hstat_elev):
        # CONVERSION
        # sun_vec az    0       -45     -90     +-180       +90     +45
        # hstat_az      -90     -45     0       +90         +-180   -135
        hstat_az = -hstat_az - N.pi/2
        
        for i in xrange(len(self.pos)):
            self.hstat_normals[i] = solar_vector(hstat_az[i], hstat_elev[i])
        print('Heliostat normals (corrected)', (N.around(self.hstat_normals, decimals=3)))
        
        global hstat_areas
        hstat_areas = [0]*len(self.pos)
        for i in xrange(len(self.pos)):
            hstat_areas[i] = (6.09**2) * abs(N.dot(-self.sun_vec, self.hstat_normals[i]))
        print('Heliostat areas seen by source', N.around(hstat_areas, decimals=4))
    
    def trace(self):
        """Generate a flux map using much more rays than drawn"""
        # Generate a large ray bundle using [a radial stagger much denser
        # than the field] a Buie sunshape with radius equal to the longest
        # dimension of the field.

        #=============
        render = False
        #=============
        
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        print('Solar vector', sun_vec)
        
        # Generate the following number of rays
	num_rays = 100000
        iters = 150

        # Initialise lists for storing receiver hit points
        xlist, ylist = [], []
        energy, pts = [], []
        
        e = TracerEngine(self.plant)

        for i in xrange(iters):
            print('ITERATION ', i+1, ' of ', iters)
                   
            # Perform the trace:
            rays = self.gen_rays(num_rays)
            a, b, hits_list = e.ray_tracer(rays, 3, 0.05, tree=True, block=False)
            del a, b
            e.minener = 1e-5

            # Render:
            if render:
                trace_scene = Renderer(e)
                trace_scene.show_rays()
        
            # Get the energy and location of all hits using optics manager
            en, pts = self.rec.get_optics_manager().get_all_hits()
            x, y = self.rec.global_to_local(pts)[:2]

            # FLUX MAP OPERATIONS
            #===========================================================================
            # Load the position and energy files for writing
            xlist = N.load('x_list.npy')
            ylist = N.load('y_list.npy')
            energy = N.load('energy_list.npy')

            # Append values from current trace
            xlist = N.append(xlist, x)
            ylist = N.append(ylist, y)
            energy = N.append(energy, en)
            print('Lengths: xlist', len(xlist), 'ylist', len(ylist))

            # Save values from current trace plus existing ones
            N.save('x_list.npy',xlist)
            N.save('y_list.npy',ylist)
            N.save('energy_list.npy',energy)

            # Reset receiver optics manager before blocking analysis
            self.rec.get_optics_manager().reset()
            #===========================================================================

            # Load the heliostat hits and theoretical hits list
            hstat_hits = N.load('hits_list.npy')
            hstat_total = N.load('total_list.npy')

            # Evaluate how many rays actually hit the heliostat in this iteration
            #===========================================================================
            # Initialise hits list for CURRENT ITERATION
            current_hits = [0] * len(self.pos)
            blocking_surfaces = []
            # Hits in stage 0
            current_hits = hits_list[0]
            print('Hits list from engine', current_hits)
            # Hits in stage 1
            print('Blocked list', hits_list[1])

            # Get number of theoretical hits on each heliostat and puts it in hstat_total
            for idx in xrange(len(self.pos)):
                # Total theoretical for the surface = area seen by source / source area * rays
                hstat_total[idx] = N.add(hstat_total[idx], hstat_areas[idx]/sun_area * num_rays)
            hstat_hits = N.add(hstat_hits, current_hits)

            blocked_vertices = N.empty((3,0))
            blocked_directions = N.empty((3,0))

            # BLOCKING OPERATIONS (iteration-wise)
            #===========================================================================
            # 1. Get indices of surfaces which have blocked some rays
            for i in xrange(len(hits_list[1])):
                if hits_list[1][i] != 0:
                    blocking_surfaces.append(i)
            print('blocking surfaces', blocking_surfaces)

            # 2. Call the optics managers of each of those surfaces and combine energies,
            #    vertices and directions into three arrays
            
            if len(blocking_surfaces) != 0:
	        for idx in blocking_surfaces:
                    energy_b, pts_b, dir_b = self.field.get_surfaces()[idx].get_optics_manager().get_all_hits()
                    blocked_vertices = N.hstack((blocked_vertices, N.array(pts_b)))
                    blocked_directions = N.hstack((blocked_directions, N.array(dir_b)))
                
                    blocked_vertices = N.array(blocked_vertices)
                    blocked_directions = N.array(blocked_directions)
                blocked_energies = [[20]] * len(blocked_vertices.T)
                blocked_energies = N.array(blocked_energies)
                del energy_b, pts_b, dir_b

            	# 3. Negate the directions array
            	blocked_directions *= -1

            	# 4. Create a new raybundle with those energues, vertices and directions
            	print('vertices', blocked_vertices)
            	print('directions', blocked_directions)
            	print('energies', blocked_energies)
            	blocked_rays = RayBundle(vertices=blocked_vertices, directions=blocked_directions, \
                	                     energy=blocked_energies)

            	# 5. Ray-trace this bundle through the scene with ONE iteration
            	a, b, blocked_list = e.ray_tracer(blocked_rays, 1, 0.05, tree=True, block=True)
            	e.minener = 1e-5
            	del a, b

            	# 6. Obtain hits_list[0] as usual and subtract this from actual hits list
            	blocked_hits = blocked_list[0]
            	hstat_hits = N.subtract(hstat_hits, blocked_hits)
            	print('>>>>>>Net hits from this iteration', hstat_hits)

            	# 7. Reset all surfaces for next iteration
            	for k in xrange(len(self.pos)):
            	    self.field.get_surfaces()[k].get_optics_manager().reset()
            #===========================================================================
        
            # Save hits info to file
            N.save('hits_list.npy', hstat_hits)
            N.save('total_list.npy', hstat_total)

        # Range of the flux map output, in metres
        rngx = 5.5
        rngy = 5.5
        
        # Flux map configuration
        bins = 50
        H, xbins, ybins = N.histogram2d(N.load('x_list.npy'), N.load('y_list.npy'), bins, \
            range=([-rngx,rngx], [-rngy,rngy]), weights=N.load('energy_list.npy'))
        extent = [ybins[0], ybins[-1], xbins[-1], xbins[0]]
        plt.imshow(H, extent=extent, interpolation='nearest')
        plt.colorbar()
        plt.title('Iteration')
        plt.show()
       
        print('Heliostat hits:', hstat_hits, 'length', len(hstat_hits))
        print('Heliostat theoretical hits:', N.ceil(hstat_total), 'length', len(hstat_total))
        print('Heliostat efficiencies', N.around((hstat_hits / N.ceil(hstat_total)), decimals=4))
        

scene = TowerScene()
hstat_az, hstat_elev = scene.aim_field()
scene.calculate_area(hstat_az, hstat_elev)
scene.trace()

