import numpy as N
from scipy.constants import degree

from tracer.models.heliostat_field import *
from tracer.models.sunshape_pillbox import *
from tracer.tracer_engine import TracerEngine
from tracer.CoIn_rendering.rendering import *


class RayTracing:
    '''
    use the TowerScene to do ray-tracing simulations
    '''

    def __init__(self, num_rays, rendering):
        '''
        num_rays - int, number of rays
        rendering - boolean, visualise the scene
        '''
        self.num_rays=num_rays
        self.rendering=rendering        
        self.gen_plant()

    def gen_plant(self):

        #--------------
        #    field
        #--------------
        # heliostats
        self.hst_w=1.85
        self.hst_h=2.44
        reflectivity=0.9
        slope_type='normal_sphere' #'normal_sphere' or 'pillbox_sphere' or 'perfect'
        slope_error=1.64e-3 #rad (0 is a perfect mirror)
        curved=True # or False: flat mirror
        self.oneMirror=False # or True: for simulating just one mirror
        hst_file='./examples/heliostat_field_example/hst_info.csv' #or None
        #hst_file=None


        if self.oneMirror:
            if hst_file!=None:
                # index of the position of the heliostats
                index=2 
                pos=N.zeros(3)
                foc=0. 
            else:
                # or hst_file=None,by putting the specific pos and focal
                index=-1 #
                pos=N.r_[0., 20., 0.] 
                foc=30.   

        else:
            index=-1
            pos=N.zeros(3)
            foc=0. 
  
        heliostat=HeliostatGenerator(self.hst_w, self.hst_h, absorptivity=1.-reflectivity, sigma_xy=slope_error, slope=slope_type,curved=curved, one_mirror=self.oneMirror, index=index, pos=pos, foc=foc)

        # layout and field

        layout=KnownField(hst_file,pos,foc)

        # tracking
        tracking_mode='TiltRoll' # 'AzEl'or 'TiltRoll'

        # aiming
        aiming_mode='SinglePoint' # 'MultiFixed' or 'SinglePoint'

        #---------------
        #   receiver
        #---------------
        rec_w=1.3
        rec_h=1.3
        absorptivity=0.96
        rec_loc=N.r_[0.,0.,26.8] # receiver location
        rec_rot=N.r_[106.,0.,0.]*degree # receiver rot through x, y, z (rad)

        receiver=FlatOneSidedReceiver(rec_w,rec_h,absorptivity)
        mount_rec=MountReceiver(rec_loc,rec_rot)

        #--------------
        #     solar
        #--------------
        sunshape='pillbox'
        self.sigma=4.65e-3
        self.DNI=1000.
        sun_az=180.
        sun_zenith=0.
        self.sun_vec=solar_vector(sun_az, sun_zenith)


        tower_scene=TowerScene(self.sun_vec, self.oneMirror)
        tower_scene(heliostat, layout, aiming_mode, tracking_mode,receiver,mount_rec)

        self.system=tower_scene.system
        self.pos=heliostat.pos


    def gen_rays(self):
        #===================================================
        if self.oneMirror:

            field_centre=self.pos
            x=self.hst_w*1.2
            y=self.hst_h*1.2
            centre = N.c_[10.*self.sun_vec + field_centre]

        else:

            t_pos = self.pos.T        
            xc_min = t_pos[0][N.argmin(t_pos[0])]
            xc_max = t_pos[0][N.argmax(t_pos[0])]
            yc_min = t_pos[1][N.argmin(t_pos[1])]
            yc_max = t_pos[1][N.argmax(t_pos[1])]

            if yc_min==yc_max:
	            y_dist=self.hst_height
            else:
	            y_dist = yc_max - yc_min
            if xc_min==xc_max:
	            x_dist=self.hst_width
            else:
	            x_dist = xc_max - xc_min

            xc_cent = (xc_min + xc_max) / 2.
            yc_cent = (yc_min + yc_max) / 2.
            field_centre = N.r_[xc_cent, yc_cent, 0.]

            # disc source covering entire field area:
            radius = 1.10 * math.sqrt((x_dist/2.)**2 + (y_dist/2.)**2)
            x=x_dist*1.2
            y=y_dist*1.2
            centre = N.c_[200.*self.sun_vec + field_centre]

        direction = N.array(-self.sun_vec)	

        rays=pillbox_rect_bundle(self.num_rays, centre, direction, x,y, self.sigma, self.DNI)
        return rays


    def trace(self):
        # define the tracer engine
        e = TracerEngine(self.system)
        e.minerer=1e-10

        #define the rays
        rays=self.gen_rays()

        # ray-tracing
        e.ray_tracer(rays, reps=100, min_energy=1e-10)

        if self.rendering:
            trace_scene=Renderer(e)
            trace_scene.show_rays(resolution=10)		     

		
if __name__=='__main__':
	
	rt=RayTracing(num_rays=10000, rendering=True)
	rt.trace()	
