'''
A tower/heliostat field example based on the Tower_gui.py example.
It has been cleaned from GUI and mayaVi stuff and given the CoIn rendering import and capacity.
'''
import traits.api as t_api
import traitsui.api as tui

from tracer.CoIn_rendering.rendering import *

import numpy as N
from scipy.constants import degree
import math

from tracer.ray_bundle import RayBundle
from tracer.sources import pillbox_sunshape_directions, buie_sunshape
from tracer.assembly import Assembly
from tracer.spatial_geometry import rotx, roty, rotation_to_z
from tracer.tracer_engine import TracerEngine

from tracer.models.one_sided_mirror import one_sided_receiver
from tracer.models.heliostat_field import HeliostatField, radial_stagger, solar_vector

class TowerScene():
    # Location of the sun:
    sun_az = 80. #80, zero is due south
    sun_elev = 45. #45, zero is up, 90 is horizon
    
    def __init__(self):
        self.gen_plant() 
   
    def gen_rays(self, num_rays):
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        centre = N.c_[20*sun_vec] # centre of the source, as a column vector, originally [0,0,10]
        x = 1/(math.sqrt(2)) # got this from /examples/test_case.py, originally 1/(math.sqrt(2))
        direction = N.array(-sun_vec) # originally [0,x,-x]

        # radius=22, CSR=2.25%, flux=1000W/m2
        rays = buie_sunshape(num_rays, centre, direction, 22, 0.0225, 1000)
        
        return rays
    
    def gen_plant(self):
        # import custom coordinate file
        self.pos = N.loadtxt("sandia_hstat_coordinates.csv", delimiter=',')
        self.pos *= 0.1
        # set heliostat field characteristics: 6.09m*6.09m, abs = 0, aim_h = 61
        self.field = HeliostatField(self.pos, 6.09e-1, 6.09e-1, 0, 6.1)

        self.rec, recobj = one_sided_receiver(1.1, 1.1) # was (1, 1)
        rec_trans = rotx(N.pi/2)
        rec_trans[2,3] = 6.1 # height of the tower
        recobj.set_transform(rec_trans)

        self.plant = Assembly(objects=[recobj], subassemblies=[self.field])
    
    def aim_field(self):
        self.field.aim_to_sun(self.sun_az*degree, self.sun_elev*degree)
    
    def trace(self):
        """Generate a flux map using much more rays than drawn"""
        # Generate a large ray bundle using [a radial stagger much denser
        # than the field] a Buie sunshape with radius equal to the longest
        # dimension of the field.
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        print(sun_vec)
        
        hstat_rays = 20
        num_rays = hstat_rays*len(self.field.get_heliostats())

        #
        #rot_sun = rotation_to_z(-sun_vec)
        #direct = N.dot(rot_sun, pillbox_sunshape_directions(num_rays, 0.00465))
        
        #xy = N.random.uniform(low=-0.25, high=0.25, size=(2, num_rays))
        #base_pos = N.tile(self.pos, (hstat_rays, 1)).T
        #base_pos += N.dot(rot_sun[:,:2], xy)
        
        #base_pos -= direct
        #
        
        #rays = RayBundle(base_pos, direct, energy=N.ones(num_rays))
        #rays = buie_sunshape(num_rays, (self.pos+sun_vec).T, N.tile(-sun_vec, (self.pos.shape[0], 1)).T, 4.65, 0.0225, 1000)
        rays = self.gen_rays(num_rays)
        
        # Perform the trace:
        e = TracerEngine(self.plant)
        e.ray_tracer(rays, 100, 0.05, tree=True)
        e.minener = 1e-5

	# Render:
        trace_scene = Renderer(e)
        trace_scene.show_rays()

scene = TowerScene()
scene.aim_field()
scene.trace()
