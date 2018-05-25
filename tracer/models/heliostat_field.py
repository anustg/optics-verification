"""
heliostat_field package development 
Since Oct 2016 

Features:
----------------------------------
20161008 -- add model of receivers
20170110 -- new receiver structure (tubes) 
20170130 -- aiming strategy
20170221 -- single mirror model
20170325 -- perfect mirror model
20170403 -- differet slope error options (pillbox, normal)
20170521 -- saving ray-bunles
20170904 -- revise AzEl Tracking
"""

from numpy import r_
import numpy as N
import types
import math 
#from ..CoIn_rendering.rendering_blades import *

from ..assembly import Assembly
from ..spatial_geometry import rotx, roty, rotz
from ..cylinder import *
from ..optics_callables import *

from ..object import AssembledObject
from ..surface import Surface
from ..flat_surface import RectPlateGM
from .. import optics_callables as opt

from .one_sided_mirror import rect_one_sided_mirror, rect_para_one_sided_mirror, perfect_rect_para_one_sided_mirror

#from ..package_dev.bladed_receiver import bladed_receiver
#from ..package_dev.bladed_receiver_tube import bladed_receiver_tube
#from ..package_dev.test_surface_reflection import receiver_test_refl


class TowerScene(Assembly):
    def __init__(self,sun_vec, one_mirror=False):
        self.one_mirror=one_mirror
        self.sun_vec=sun_vec

    def __call__(self, field, layout, aiming_mode, tracking_mode,receiver, mount_rec):

        if field == None:
            # use the saved ray bundles 
            mount_rec(receiver._rec)								
            self.system=Assembly(objects=[receiver._rec])	

        else:
            field(layout)
            heliostats=Assembly(objects=field._heliostats)

            # Aiming
            if aiming_mode == 'MultiFixed':
                aiming=MultiFixedPointsAiming(layout.aiming)
            elif aiming_mode =='SinglePoint':
                aiming=SinglePointAiming(layout.pos, mount_rec.position, self.one_mirror)

            # Tracking
            if tracking_mode =='TiltRoll':
                tracking=TiltRoll(self.sun_vec,self.one_mirror)	
            elif tracking_mode =='AzEl':
                tracking=AzElTrackings(self.sun_vec,self.one_mirror) 

            if self.one_mirror:
                tracking.aim_to_sun(-field.pos, aiming.aiming_points, self.one_mirror)
            else:		
                tracking.aim_to_sun(-layout.pos, aiming.aiming_points, self.one_mirror)
            tracking(field)

            if receiver==None:
                # for saving ray bundles
                self.system=heliostats
            else:
                mount_rec(receiver._rec)
                self.system=Assembly(objects=[receiver._rec], subassemblies=[heliostats])


class MountReceiver:		
    def __init__(self,position, rotation):

        '''
        Arguments:
        position: receiver position (3,) array
	        position[0]-x
	        position[1]-y
	        position[2]-z
        rotation: (3,) array, contain the rotation angle (randian) 
	        rotation[0]-x_ang
	        rotation[1]-y_ang
	        rotation[2]-z_ang
        '''
        self.position=position
        self.rotation=rotation

    def __call__(self, receiver):

        trans =N.dot(rotx(self.rotation[0]),rotz(self.rotation[2]))
        trans[0,3]=self.position[0]
        trans[1,3]=self.position[1]
        trans[2,3]=self.position[2]
        receiver.set_transform(trans)

       
class Receiver:
    def __init__(self):
        pass

class FlatPlate(Receiver):
    def __init__(self, width, height,absorptivity):
        """
        construct a flat one surface on the XY plane, 
        which is a ReflctiveReceiver object

        Arguments:
	        width - the extent along the x axis in the local frame.
	        height - the extent along the y axis in the local frame.
	        absorptivity - the ratio of energy incident on the reflective side that's
	         not reflected back.

        Returns:
        obj-The AssembledObject	

        """
        surface=Surface(RectPlateGM(width, height), 
		        opt.ReflectiveReceiver(absorptivity))

        obj = AssembledObject(surfs=[surface])

        #obj.surfaces_for_next_iteration = types.MethodType(
            #surfaces_for_next_iteration, obj, obj.__class__)

        self._rec=obj

class FlatOneSidedReceiver(Receiver):
    def __init__(self, width, height,absorptivity):

        """
	        construct an object with two surfaces: one on the XY plane, that is
	        specularly reflective, and one slightly below (negative z), that is opaque.
	        The reflective surface is by default also apaque; it is a ReflectiveReceiver
	        object, so all hits can be obtained after a trace using that surface's 
	        get_all_hits() method.

        Arguments:
	        width - the extent along the x axis in the local frame.
	        height - the extent along the y axis in the local frame.
	        absorptivity - the ratio of energy incident on the reflective side that's
	         not reflected back.

         Returns:
	        front - the receiving surface
	        obj - the AssembledObject containing both surfaces
         """

        front = Surface(RectPlateGM(width, height), 
		        opt.ReflectiveReceiver(absorptivity))
        back = Surface(RectPlateGM(width, height), opt.ReflectiveReceiver(1.),
		        location=r_[0., 0., -1e-10])
        obj = AssembledObject(surfs=[front, back])

        #obj.surfaces_for_next_iteration = types.MethodType(
            #surfaces_for_next_iteration, obj, obj.__class__)

        self._rec=obj


class BladedReceiver(Receiver):
    def __init__(self, W, H, D, T, a, num, edge1, edge2,absorptivity,bld_ends, spacespecial):
        '''
        Build the bladed receiver

        Arguments:
        W-- width of blade along the x axis, m    
        H-- height of the bladed receiver along the z axis, m
        D-- depth of blade along the y axis, m
        T-- wall thickness, m
        a-- angle between blades and the back wall, radian 
        num-- number of blades
        edge1 -- width of the edge panelAB
        edge2 -- height of the edge panelCD
        absorptivity-absorptivity of the blades surfaces
        bld_ends -- bool, whether to have a side wall to block the blade ends
        spacespecial -- bool, whether has the extend space to make (num of space = num of blade)

        Returns:
        obj -The AssembledObject
        '''  
        obj =AssembledObject(surfs=bladed_receiver(W, H, D, T, a, num, edge1, edge2, absorptivity,bld_ends,spacespecial))

        #obj.surfaces_for_next_iteration = types.MethodType(
            #surfaces_for_next_iteration, obj, obj.__class__)

        self._rec=obj

class CylinderReceiver(Receiver):
    def __init__(self, h, d, absorptivity):
        '''	
        build a big cylinder to test how to implement tubes for bladed receiver
        h: cylinder height z (m)
        d: cylinder diameter (m)  
        '''
        obj= AssembledObject(surfs=[Surface(FiniteCylinder(diameter=d, height=h),LambertianReceiver(absorptivity))])

        self._rec=obj

class BladedTubeReceiver(Receiver):
    def __init__(self,h, d, nt_bkw, nt_bl, nt_edge, num, absorptivity):
        '''
        bladed receiver with tubes
        (now the blades perpendicular to back wall)

        Arguments:
        h -- the height of the backwall along the z axis, m
        d -- tube diameter, m
        nt_bkw -- number of tubes per space on the back wall
        nt_bl -- number of tubes per blade
        nt_edge -- number of tubes for the edge section (now edge1=edge2)
        num-- number of blades
        absorptivity-absorptivity of the blades surfaces

        '''

        obj =AssembledObject(surfs=bladed_receiver_tube(h, d, nt_bkw, nt_bl, nt_edge, num, absorptivity))

        self._rec=obj	

class FlatTubeReceiver(Receiver):
    def __init__(self, h, d, nt, absorptivity):
        '''
        flat receiver with tubes

        Arguments:
        h: tube height, the height of the flat receiver along the z axis, m
        d: tube diameter, m
        nt: number of tubes
        absorptivity -absorptivity of the tube surfaces
        '''	
        H=h
        W= nt*d 
        surfs=[]
        for i in xrange(int(nt)):
            tube=Surface(FiniteCylinder(diameter=d, height=h),LambertianReceiver(absorptivity))	
            trans=N.eye(4)
            trans[0,3]=W/2.-d/2.-i*d
            tube.set_transform(trans)
            surfs.append(tube)
        obj=AssembledObject(surfs=surfs)

        self._rec=obj


class ReceiverTestRefl(Receiver):
    def __init__(self,absorptivity):
        '''
        Build two flat surfaces which can reflect rays to each other easily

        Arguments:

        absorptivity-absorptivity of the surfaces
        '''     
        obj =AssembledObject(surfs=receiver_test_refl(absorptivity))
        #obj.surfaces_for_next_iteration = types.MethodType(
            #surfaces_for_next_iteration, obj, obj.__class__)

        self._rec=obj	

class FieldLayout:
    def __init__(self, positions,focal_lengths=None):
        '''
        Creat heliostat field layout, 
        including the position and focal length of each heliostat.

        Arguments:
        positions - an (n,3) array, each row has the location of one heliostat
        focal_lengthes - an (n,1) array, each element is the focal length of one heliostat
        '''

        self.pos=positions
        self.focals=focal_lengths

    def get_positions(self):
        return self.pos

    def get_focals(self):
        return self.focals

    def get_one_mirror(self,index):
        '''
        Arguments:
        index - index of the selected mirror in the field
        '''
        one_pos=self.pos[index]
        one_foc=self.focals[index]

        return one_pos, one_foc


class KnownField(FieldLayout):
    def __init__(self, filename,pos=N.zeros(3),foc=N.zeros(3)):
        '''
        For already known field layouts.

        Arguments:
        filename - filename of the file (.csv) that listed the coordinates and 
            focal lengths of the known field
	            layout[:,0]- x coordinates
	            layout[:,1]- y coordinates
	            layout[:,2]- z coordinates
	            layout[:,3]- focal lengthes 			

        self.pos -An array with an x,y,z row for each heliostat (shape n,3)
        self.focals- An array with focal length for each heliostat (shape n,1)
        '''
        if (any(pos) == 0. and any(foc)==0.):
            info=N.loadtxt(filename, delimiter=',', skiprows=1)
            self.pos=info[:,:3]
            self.focals=info[:,3]

            if len(info[0,:])>4:
                self.aiming=info[:,4:]

        else:
            self.pos=pos
            self.focals=foc
				
		
class RadialStaggerLT(FieldLayout):

    def __init__(self, start_ang, end_ang, az_space, rmin, rmax, r_space,z, focal_lengths=None):
        """
        Calculate positions of heliostats in a radial-stagger field. This is a
        common way to arrange heliostats.

        Arguments:
        start_ang, end_ang - the angle in radians CW from the X axis that define
        the field's boundaries.
        az_space - the azimuthal space between two heliostats, in [rad]
        rmin, rmax - the boundaries of the field in the radial direction.
        r_space - the space between radial lines of heliostats.
        z - the height of the heliostat to the ground

        Returns:
        An array with an x,y,z row for each heliostat (shape n,3)
        """		
        self.start_ang=start_ang
        self.end_ang=end_ang
        self.za_space=az_space
        self.rmin=rmin
        self.rmax=rmax
        self.r_space=r_space
        self.z=z

        rs = N.r_[self.rmin:self.rmax:self.r_space]
        angs = N.r_[self.start_ang:self.end_ang:self.az_space/2]

        # 1st stagger:
        xs1 = N.outer(rs[::2], N.cos(angs[::2])).flatten()
        ys1 = N.outer(rs[::2], N.sin(angs[::2])).flatten()

        # 2nd staggeer:
        xs2 = N.outer(rs[1::2], N.cos(angs[1::2])).flatten()
        ys2 = N.outer(rs[1::2], N.sin(angs[1::2])).flatten()

        xs = N.r_[xs1, xs2]
        ys = N.r_[ys1, ys2]

        zs=N.ones(N.shape(xs))*self.z

        self.pos=N.vstack((xs, ys,zs)).T

        self.focals=focal_lengths


class HeliostatGenerator:

    def __init__(self, width, height, absorptivity, sigma_xy,slope, curved=True, one_mirror=False, index=-1,pos=N.zeros(3),foc=N.zeros(3)):
        """
        Generates heliostats, each being a rectangular one-sided
        mirror (either curved or flat), initially pointing downward 
        - for safety reasons, of course :)

        Arguments:		
        width, height - The width and height, respectively, of each heliostat.
        absorptivity - part of incident energy absorbed by the heliostat.
        sigma_xy: slope error (mrad)
        slope: the shape of slope error, options: (1)pillbox (2) normal (3) ns (normal simplified)
        curved: boolean, True or False  
        if simulate one mirror:
        index: index of the selected mirror in the field
        pos, foc: position and focal length of the selected mirror
        """

        self.width=width
        self.height=height
        self.absorpt=absorptivity
        self.sigma=sigma_xy
        self.slope=slope
        self.curved=curved
        self.one_mirror=one_mirror
        if self.one_mirror:
            self.index=index
            self.pos=pos
            self.foc=foc


    def __call__(self, layout):
		
        face_down = rotx(N.pi)
        self._heliostats = []

        if self.one_mirror:
            if self.index !=-1:
                self.pos, self.foc=layout.get_one_mirror(self.index)	
        else:
            self.pos=layout.pos


        if self.curved:
            if self.one_mirror:
                if self.sigma==0. or self.slope=='perfect':
                    hstat = perfect_rect_para_one_sided_mirror(self.width, self.height, self.foc, self.absorpt)
                else:
                    hstat = rect_para_one_sided_mirror(self.width, self.height, self.foc, self.absorpt, self.sigma,self.slope)				  
                trans = face_down.copy()
                trans[:3,3] = self.pos
                hstat.set_transform(trans)
                self._heliostats.append(hstat)			

            else:
                print "I'm generating curved mirrors"
                for i in xrange(len(self.pos[:,0])):
                    hstat_pos = N.array(self.pos[i])			
                    focal= layout.focals[i]
                    hstat = rect_para_one_sided_mirror(self.width, self.height, focal, self.absorpt, self.sigma,self.slope)	

                    trans = face_down.copy()
                    trans[:3,3] = hstat_pos 
                    hstat.set_transform(trans)
                    self._heliostats.append(hstat)

        else:
            if self.one_mirror:
                hstat =rect_one_sided_mirror(self.width, self.height,self.absorpt, self.sigma,self.slope)
                trans = face_down.copy()
                trans[:3,3] = self.pos 
                hstat.set_transform(trans)
                self._heliostats.append(hstat)	

            else:
                for pos in self.pos:
                    hstat_pos = N.array(pos)			
                    hstat =rect_one_sided_mirror(self.width, self.height,self.absorpt, self.sigma,self.slope)	
                    trans = face_down.copy()
                    trans[:3,3] = pos
                    hstat.set_transform(trans)
                    self._heliostats.append(hstat)

    def get_heliostats(self):
        '''
        Access the list of the mirrors representing the heliostats
        '''
        return self._heliostats

    def get_one_mirror(self):
        return self.pos, self.foc


class AimingStrategy():
    def __init__(self):

        '''
        Define aiming strategy: single aiming or multiple aiming points

        return the array of aiming points- - an (n,3) array, each row represent the aiming point for the corresponding heliostat 
        '''

        pass

class SinglePointAiming(AimingStrategy):

    def __init__(self, hstat_pos, single_point,one_mirror=False):

        '''
        Arguement:
        hstat_pos- an (n,3) array, each row has the location of one heliostat
        single_point - an (1,3) array, represent the single aiming point

        Return:
        aiming_points - an (n,3) array, each row represent the aiming point for 
	        the corresponding heliostat 
        '''

        if one_mirror:
            self.aiming_points=single_point

        else:
            self.aiming_points=N.zeros(N.shape(hstat_pos))

            self.aiming_points[:,0]=single_point[0]		
            self.aiming_points[:,1]=single_point[1]
            self.aiming_points[:,2]=single_point[2]


class MultiFixedPointsAiming(AimingStrategy):

    def __init__(self,aiming_points):
        '''
        For already known field layouts and aiming points.

        Arguments:
        aiming_points - an (n,3) array, each row represent the aiming point
	        for the corresponding helisotat
        '''
        self.aiming_points=aiming_points


class SimpleLinkedPointsAiming(AimingStrategy):

    def __init__(self, hstat_pos, tgt_pos, tgt_width, tgt_height,rot_angle,ap_x,ap_z,spil_contrl):
        '''
        Simply assign each individul mirror to an aiming point on the target
        eg: if the field is a rectancular, the traget aperture is a square
        farther mirrors aim to the top of the square aperture (left to left)
        nearer mirrors aim to the bottom of the square aperture (right to right)


        Arguments:
        hstat_pos- an (n,3) array, each row has the location of one heliostat		
        tgt_pos-  an (1,3) array, represent the target position
        tgt_width- the width of the aperture of the target (x direction)
        tgt_height- the height of the aperture of the target (z direction)
        rot_angle- target rotating angle, in radius: 
					            pi/2   -> vertical
					            pi*3/4 -> 45 dg face to the field
					            26+90dg-> 26 dg between the z axies
        ap_x-aiming points in the x directions 
        ap_z-aiming points in the z directions
        spil_contrl- to control the spillage, distance to the edge of the traget (m)

        self.aiming_points - an (n,3) array, each row represent the aiming point
	        for the corresponding helisotat
        '''
        # the heliosts at edges
        x_max=N.max(hstat_pos[:,0])
        x_min=N.min(hstat_pos[:,0])
        y_max=N.max(hstat_pos[:,1])
        y_min=N.min(hstat_pos[:,1])


        # rotating angle
        theta=rot_angle-N.pi/2.

        # edges of the aiming points on the traget
        xt_max=tgt_width/2.-spil_contrl
        xt_min=-tgt_width/2.+spil_contrl
        zt_max=(tgt_height/2.-spil_contrl)*N.cos(theta)+tgt_pos[2]
        zt_min=(-tgt_height/2.+spil_contrl)*N.cos(theta)+tgt_pos[2]


        # aiming points on the traget
        distance_x=(xt_max-xt_min)/(ap_x-1)
        distance_z=(zt_max-zt_min)/(ap_z-1)
        aim_x=N.linspace(xt_min,xt_max+distance_x,ap_x)
        aim_z=N.linspace(zt_min,zt_max+distance_z,ap_z)


        # grids for the field
        gd_x=N.linspace(x_min-5.,x_max+5.,ap_x)
        gd_y=N.linspace(y_min-5.,y_max+5.,ap_z)



        self.aiming_points=N.zeros(N.shape(hstat_pos))

        for h in xrange(len(hstat_pos[:,0])):

            for i in xrange(ap_x-1):
                for j in xrange(ap_z-1):
                    if hstat_pos[h,0]>=gd_x[i] and hstat_pos[h,0]<gd_x[i+1] and hstat_pos[h,1]>=gd_y[j] and hstat_pos[h,1]<gd_y[j+1]:
                        self.aiming_points[h,0]=aim_x[i+1]
                        self.aiming_points[h,2]=aim_z[-j-1]
                        self.aiming_points[h,1]=N.tan(theta)*(aim_z[-j-1]-tgt_pos[2])
						

class TrackingSystem:
    def __init__(self,sun_vec):		
        '''
        implement two different tracking system, namely:
        (1) AzEl- Azimuth-elevation tracking
        (2) TiltRoll -pitch/tilt-roll tracking

        Argument:
        sun_vec : a 3-component 1D array with the solar vector.
        '''	
        self.sun_vec=sun_vec

    def aim_to_sun(self, tower_vec, aiming_pos, one_mirror=False):
        """
        Aim the heliostats in a direction that brings the incident energy to
        the aiming points.
	
        return: hstat_norm - an (n,3) array, each row has the normal vector of 
	        each heliostat when it aims to sun

        """
        if one_mirror:
            tower_vec += aiming_pos
            tower_vec /= N.sqrt(N.sum(tower_vec**2))
            self.hstat_norm=self.sun_vec+tower_vec
            self.hstat_norm /= N.sqrt(N.sum(self.hstat_norm**2))


        else:
            tower_vec += aiming_pos
            tower_vec /= N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])
            self.hstat_norm=self.sun_vec+tower_vec
            self.hstat_norm /= N.sqrt(N.sum(self.hstat_norm**2, axis=1)[:,None])


class AzElTrackings(TrackingSystem):
    def __init__(self,sun_vec,one_mirror=False):
        '''
        Azimuth and elevation tracking

        Argument:
        sun_vec : a 3-component 1D array with the solar vector.

        '''

        self.sun_vec=sun_vec
        self.one_mirror=one_mirror

    def __call__(self,heliostats):

        if self.one_mirror:

            hstat_elev = N.arccos(self.hstat_norm[2])

            norm_x=self.hstat_norm[0]
            norm_y=self.hstat_norm[1]
            norm_z=self.hstat_norm[2]

            if norm_x>=0:
                hstat_az = N.arccos(-norm_y/N.sqrt(norm_x**2+norm_y**2))                                     
            elif norm_x<0:
                hstat_az = N.arccos(norm_y/N.sqrt(norm_x**2+norm_y**2)) +N.pi
 	
            az_rot = rotz(hstat_az)
            elev_rot = rotx(hstat_elev)

            trans = N.dot(az_rot, elev_rot)
            trans[:3,3] = heliostats.pos          
            heliostats._heliostats[0].set_transform(trans)

        else:
            norm_x=self.hstat_norm[:,0]
            norm_y=self.hstat_norm[:,1]
            norm_z=self.hstat_norm[:,2]
						
            hstat_elev = N.arccos(norm_z)

            for hidx in xrange(heliostats.pos.shape[0]):
                if norm_x[hidx]>=0:
                    hstat_az = N.arccos(-norm_y[hidx]/N.sqrt(norm_x[hidx]**2+norm_y[hidx]**2))                                     
                elif norm_x[hidx]<0:
                    hstat_az = N.arccos(norm_y[hidx]/N.sqrt(norm_x[hidx]**2+norm_y[hidx]**2)) +N.pi
     	
                elev_rot = rotx(hstat_elev[hidx])
                az_rot = rotz(hstat_az)

                trans = N.dot(az_rot,elev_rot)

                trans[:3,3] = heliostats.pos[hidx]

                heliostats._heliostats[hidx].set_transform(trans)


class TiltRoll(TrackingSystem):
    def __init__(self,sun_vec,one_mirror=False):
        '''
        Pitch and roll tracking

        Argument:
        sun_vec : a 3-component 1D array with the solar vector.

        hstat_alfa - the angle that the mirror should rotat 
			        around x-axis to aim to the sun
        hstat_beta - the angle that the mirror should rotat 
			        around y-axis to aim to the sun
        '''
        self.sun_vec=sun_vec 
        self.one_mirror=one_mirror

    def __call__(self,heliostats):


        if self.one_mirror:
	
            '''
            #angle between x-axis and the normal vector of the mirror when it aims to sun
            gama=N.arccos(self.hstat_norm[0]) 

            if self.hstat_norm[1]>0:			
	            hstat_alfa=-N.arccos(1.-((N.sin(gama)-self.hstat_norm[2])**2+self.hstat_norm[1]**2)/(2.*(N.sin(gama))**2)) 
            else:
	            hstat_alfa=N.arccos(1.-((N.sin(gama)-self.hstat_norm[2])**2+self.hstat_norm[1]**2)/(2.*(N.sin(gama))**2)) 

            hstat_beta=N.arcsin(self.hstat_norm[0])
            '''

            hstat_alfa=-N.arctan2(self.hstat_norm[1],self.hstat_norm[2])
            hstat_beta=N.arcsin(self.hstat_norm[0])
            print hstat_alfa
            print self.hstat_norm[1]
            print self.hstat_norm[2]

            alfa_rot = rotx(hstat_alfa)
            beta_rot = roty(hstat_beta)

            trans = N.dot(alfa_rot,beta_rot)
            trans[:3,3] = heliostats.pos

            heliostats._heliostats[0].set_transform(trans)


        else:
            #angle between x-axis and the normal vector of the mirror when it aims to sun
            #gama=N.arccos(self.hstat_norm[:,0])  
            #hstat_alfa=-N.arccos(1.-((N.sin(gama)-self.hstat_norm[:,2])**2+self.hstat_norm[:,1]**2)/(2.*(N.sin(gama))**2)) 

            hstat_alfa=-N.arctan2(self.hstat_norm[:,1],self.hstat_norm[:,2])
            hstat_beta=N.arcsin(self.hstat_norm[:,0])

            for hidx in xrange(heliostats.pos.shape[0]):
                alfa_rot = rotx(hstat_alfa[hidx])
                beta_rot = roty(hstat_beta[hidx])

                trans = N.dot(alfa_rot, beta_rot)
                trans[:3,3] = heliostats.pos[hidx]

                heliostats._heliostats[hidx].set_transform(trans)

def solar_vector(azimuth, zenith):
    """
    Calculate the solar vector using elevation and azimuth.

    Arguments:
    azimuth - the sun's azimuth, in radians, 
	    from South increasing towards to the West
    zenith - angle created between the solar vector and the Z axis, in radians.

    Returns: a 3-component 1D array with the solar vector.
    """
    sun_z = N.cos(zenith)
    sun_y=-N.sin(zenith)*N.cos(azimuth)
    sun_x=-N.sin(zenith)*N.sin(azimuth)
    sun_vec = N.r_[sun_x, sun_y,sun_z] 

    return sun_vec


def surfaces_for_next_iteration(self, rays, surface_id):
    """
    Informs the ray tracer that some of the surfaces can be skipped in the
    next ireration for some of the rays.
    This implementation marks all surfaces as irrelevant to all rays.

    Arguments:
    rays - the RayBundle to check. 
    surface_id - the index of the surface which generated this bundle.

    Returns:
    an array of size s by r for s surfaces in this object and r rays,
        stating whether ray i=1..r should be intersected with surface j=1..s
        in the next iteration.
    """
    return N.zeros((len(self.surfaces), rays.get_num_rays()), dtype=N.bool)

