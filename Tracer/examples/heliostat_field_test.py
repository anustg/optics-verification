'''
A tower/heliostat field example based on the Tower_gui.py example.
It has been cleaned from GUI and mayaVi stuff and given the CoIn rendering import and capacity.
'''
import traits.api as t_api
import traitsui.api as tui

from tracer.CoIn_rendering.rendering import *

import numpy as N
from scipy.constants import degree

from tracer.ray_bundle import RayBundle
from tracer.sources import pillbox_sunshape_directions
from tracer.assembly import Assembly
from tracer.spatial_geometry import roty, rotx, rotation_to_z
from tracer.tracer_engine import TracerEngine

from tracer.models.one_sided_mirror import one_sided_receiver
from tracer.models.heliostat_field import HeliostatField, radial_stagger, solar_vector

class TowerScene():
    # Location of the sun:
    sun_az = 80. # degrees from positive X-axis
    sun_elev = 45. # degrees from XY-plane
    
    def __init__(self):
        self.gen_plant() 
   
    def gen_rays(self):
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        rpos = (self.pos + sun_vec).T
        direct = N.tile(-sun_vec, (self.pos.shape[0], 1)).T
        rays = RayBundle(rpos, direct, energy=N.ones(self.pos.shape[0]))
        
        return rays
    
    def gen_plant(self):
        # import custom coordinate file
        self.pos = N.loadtxt("sandia_hstat_coordinates.csv", delimiter=',')
        self.pos *= 0.1
        # set heliostat field characteristics: 0.52m*0.52m, abs = 0, aim_h = 61
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
        # Generate a large ray bundle using a radial stagger much denser
        # than the field.
        sun_vec = solar_vector(self.sun_az*degree, self.sun_elev*degree)
        
        hstat_rays = 20
        num_rays = hstat_rays*len(self.field.get_heliostats())
        rot_sun = rotation_to_z(-sun_vec)
        direct = N.dot(rot_sun, pillbox_sunshape_directions(num_rays, 0.00465))
        
        xy = N.random.uniform(low=-0.25, high=0.25, size=(2, num_rays))
        base_pos = N.tile(self.pos, (hstat_rays, 1)).T
        base_pos += N.dot(rot_sun[:,:2], xy)
        
        base_pos -= direct
        rays = RayBundle(base_pos, direct, energy=N.ones(num_rays))
        
        # Perform the trace:
        e = TracerEngine(self.plant)
        e.ray_tracer(rays, 100, 0.05, tree=True)
        e.minener = 1e-6 # default 1e-5

	# Render:
        trace_scene = Renderer(e)
        trace_scene.show_rays()

scene = TowerScene()
scene.aim_field()
scene.trace()
