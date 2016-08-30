"""
Manage a field of flat-surface heliostats aimed at a solar tower. The tower
is assumed to be at the origin, and the aiming is done by giving the sun's
azimuth and elevation.

The local coordinates system assumes that +x is north and +z is up, which is
also known as the Local Coordinates System in solar energy.

References:
.. [1] http://www.flickr.com/photos/8242576@N06/2652388885


#######################################################
(1)change the def solar_vertor() and def aim_to_sun():
elevation to zenith--Jan 2016

(2) add a new class--curved heliostat mirror--9 June 2016

by Ye

"""

import numpy as N

from ..assembly import Assembly
from .one_sided_mirror import rect_one_sided_mirror, rect_para_one_sided_mirror
from ..spatial_geometry import rotx, roty, rotz

class HeliostatField(Assembly):
    def __init__(self, positions, width, height, absorptivity, aim_height, sigma_xy):
        """
        Generates a field of heliostats, each being a rectangular one-sided
        mirror, initially pointing downward - for safety reasons, of course :)
        This setting is used in the Weizmann Institute Tower [1], among others.
        
        Arguments:
        positions - an (n,3) array, each row has the location of one heliostat.
        width, height - The width and height, respectively, of each
            heliostat.
        apsorpt - part of incident energy absorbed by the heliostat.
        aim_height - the height (Z coordinate) of the target for aiming
        """
        self._pos = positions  # Save collecting positions from the hstats.
        self._th = aim_height
        face_down = rotx(N.pi)
        tower_ht = N.array([0,0,self._th]) #
        
        self._heliostats = []
        for pos in positions:
            hstat_pos = N.array(pos)
            focal_length = N.linalg.norm(tower_ht - hstat_pos)
            #print(focal_length)
          	# print('initialising parabolic mirror here')
            hstat = rect_one_sided_mirror(width, height, absorptivity, sigma_xy)
            #hstat = rect_para_one_sided_mirror(width, height, focal_length, absorpt, sigma_xy)
            #print('finished initialising heliostat')
            trans = face_down.copy()
            trans[:3,3] = pos
            hstat.set_transform(trans)
            self._heliostats.append(hstat)
            #print('appended heliostat')
            
        Assembly.__init__(self, objects=self._heliostats)
    
    def get_heliostats(self):
        """Access the list of one-sided mirrors representing the heliostats"""
        return self._heliostats
    
    def set_aim_height(self, h):
        """Change the verical position of the tower's target."""
        self._th = h
    
    def aim_to_sun(self, azimuth, zenith):
        """
        Aim the heliostats in a direction that brings the incident energy to
        the tower.
        
        Arguments:
        azimuth - the sun's azimuth, in radians from east, counterclockwise.
        zenith - angle created between the solar vector and the Z axis, 
            in radians.
        """
        sun_vec = solar_vector(azimuth, zenith)
        tower_vec = -self._pos 
        tower_vec[:,2] += self._th
        tower_vec /= N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])
        hstat = sun_vec + tower_vec
        hstat /= N.sqrt(N.sum(hstat**2, axis=1)[:,None])
        
        hstat_az = N.arctan2(hstat[:,1], hstat[:,0])
        hstat_elev = N.arccos(hstat[:,2])
        
        for hidx in xrange(self._pos.shape[0]):
            az_rot = rotz(hstat_az[hidx])
            elev_rot = roty(hstat_elev[hidx])
            
            trans = N.dot(az_rot, elev_rot)
            trans[:3,3] = self._pos[hidx]
            
            self._heliostats[hidx].set_transform(trans)

	return hstat_az, hstat_elev

def solar_vector(azimuth, zenith):
    """
    Calculate the solar vector using elevation and azimuth.
    
    Arguments:
    azimuth - the sun's azimuth, in radians, from South increasing towards to the West
    zenith - angle created between the solar vector and the Z axis, 
        in radians.
    
    Returns: a 3-component 1D array with the solar vector.
    """
    sun_z = N.cos(zenith)
    sun_y=-N.sin(zenith)*N.cos(azimuth)
    sun_x=-N.sin(zenith)*N.sin(azimuth)
    sun_vec = N.r_[sun_x, sun_y,sun_z] 

    return sun_vec

def radial_stagger(start_ang, end_ang, az_space, rmin, rmax, r_space):
    """
    Calculate positions of heliostats in a radial-stagger field. This is a
    common way to arrange heliostats.
    
    Arguments:
    start_ang, end_ang - the angle in radians CW from the X axis that define
        the field's boundaries.
    az_space - the azimuthal space between two heliostats, in [rad]
    rmin, rmax - the boundaries of the field in the radial direction.
    r_space - the space between radial lines of heliostats.
    
    Returns:
    An array with an x,y row for each heliostat (shape n,2)
    """
    rs = N.r_[rmin:rmax:r_space]
    angs = N.r_[start_ang:end_ang:az_space/2]
    
    # 1st stagger:
    xs1 = N.outer(rs[::2], N.cos(angs[::2])).flatten()
    ys1 = N.outer(rs[::2], N.sin(angs[::2])).flatten()
    
    # 2nd staggeer:
    xs2 = N.outer(rs[1::2], N.cos(angs[1::2])).flatten()
    ys2 = N.outer(rs[1::2], N.sin(angs[1::2])).flatten()
    
    xs = N.r_[xs1, xs2]
    ys = N.r_[ys1, ys2]
    
    return N.vstack((xs, ys)).T


class HeliostatCurved(Assembly):
    def __init__(self, positions, width, height, focal_length_mirr, absorptivity, aim_height, sigma_xy):
        """
        Generates a field of heliostats, each being a rectangular one-sided
        mirror, initially pointing downward - for safety reasons, of course :)
        This setting is used in the Weizmann Institute Tower [1], among others.
        
        Arguments:
        positions - an (n,3) array, each row has the location of one heliostat.
        width, height - The width and height, respectively, of each
            heliostat.
        apsorpt - part of incident energy absorbed by the heliostat.
        aim_height - the height (Z coordinate) of the target for aiming
        """
        self._pos = positions  # Save collecting positions from the hstats.
        self._th = aim_height
        face_down = rotx(N.pi)
        tower_ht = N.array([0.,0.,self._th]) #
        
        self._heliostats = []
        for pos in positions:
            hstat_pos = N.array(pos)
            focal_length= N.linalg.norm(tower_ht - hstat_pos)
            #print(focal_length)
          	# print('initialising parabolic mirror here')
            hstat = rect_para_one_sided_mirror(width, height, focal_length_mirr, absorptivity, sigma_xy)
            #hstat = rect_para_one_sided_mirror(width, height, focal_length, absorpt, sigma_xy)
            #print('finished initialising heliostat')
            trans = face_down.copy()
            trans[:3,3] = pos
            hstat.set_transform(trans)
            self._heliostats.append(hstat)
            #print('appended heliostat')
            
        Assembly.__init__(self, objects=self._heliostats)
    
    def get_heliostats(self):
        """Access the list of one-sided mirrors representing the heliostats"""
        return self._heliostats
    
    def set_aim_height(self, h):
        """Change the verical position of the tower's target."""
        self._th = h
    
    def aim_to_sun(self, azimuth, zenith):
        """
        Aim the heliostats in a direction that brings the incident energy to
        the tower.
        
        Arguments:
        azimuth - the sun's azimuth, in radians from east, counterclockwise.
        zenith- angle created between the solar vector and the Z axis, 
            in radians.
        """
        sun_vec = solar_vector(azimuth, zenith)
        tower_vec = -self._pos 
        tower_vec[:,2] += self._th
        tower_vec /= N.sqrt(N.sum(tower_vec**2, axis=1)[:,None])
        hstat = sun_vec + tower_vec
        hstat /= N.sqrt(N.sum(hstat**2, axis=1)[:,None])
        
        hstat_az = N.arctan2(hstat[:,1], hstat[:,0])
        hstat_elev = N.arccos(hstat[:,2])
        
        for hidx in xrange(self._pos.shape[0]):
            az_rot = rotz(hstat_az[hidx])
            elev_rot = roty(hstat_elev[hidx])
            
            trans = N.dot(az_rot, elev_rot)
            trans[:3,3] = self._pos[hidx]
            
            self._heliostats[hidx].set_transform(trans)

	return hstat_az, hstat_elev


class one_HeliostatCurved(Assembly):
    def __init__(self, positions, width, height, focal_length_mirr, absorptivity, aim_height, sigma_xy):
        """
        Generates a field of heliostats, each being a rectangular one-sided
        mirror, initially pointing downward - for safety reasons, of course :)
        This setting is used in the Weizmann Institute Tower [1], among others.
        
        Arguments:
        positions - an (n,3) array, each row has the location of one heliostat.
        width, height - The width and height, respectively, of each
            heliostat.
        apsorpt - part of incident energy absorbed by the heliostat.
        aim_height - the height (Z coordinate) of the target for aiming
        """
        self._pos = positions  # Save collecting positions from the hstats.
        self._th = aim_height
        face_down = rotx(N.pi)
        tower_ht = N.array([0.,0.,self._th]) #
        
        self._heliostats = []
        for pos in positions:
            hstat_pos = N.array(pos)
            focal_length= N.linalg.norm(tower_ht - hstat_pos)
            #print(focal_length)
          	# print('initialising parabolic mirror here')
            hstat = rect_para_one_sided_mirror(width, height, focal_length_mirr, absorptivity, sigma_xy)
            #hstat = rect_para_one_sided_mirror(width, height, focal_length, absorpt, sigma_xy)
            #print('finished initialising heliostat')
            trans = face_down.copy()
            trans[:3,3] = pos
            hstat.set_transform(trans)
            self._heliostats.append(hstat)
            #print('appended heliostat')
            
        Assembly.__init__(self, objects=self._heliostats)
    
    def get_heliostats(self):
        """Access the list of one-sided mirrors representing the heliostats"""
        return self._heliostats
    
    def set_aim_height(self, h):
        """Change the verical position of the tower's target."""
        self._th = h
    
    def aim_to_sun(self, azimuth, zenith):
        """
        Aim the heliostats in a direction that brings the incident energy to
        the tower.
        
        Arguments:
        azimuth - the sun's azimuth, in radians from east, counterclockwise.
        zenith - angle created between the solar vector and the Z axis, 
            in radians.
        """
        sun_vec = solar_vector(azimuth, zenith)
        tower_vec = -self._pos 
        tower_vec[2] += self._th
        tower_vec /= N.linalg.norm(tower_vec)
        hstat = sun_vec + tower_vec
        hstat /= N.linalg.norm(hstat)
        
        hstat_az = N.arctan2(hstat[1], hstat[0])
        hstat_elev = N.arccos(hstat[2])
        
	az_rot = rotz(hstat_az)
	elev_rot = roty(hstat_elev)
            
	trans = N.dot(az_rot, elev_rot)
	trans[:3,3] = self._pos
            
	self._heliostats[0].set_transform(trans)

	return hstat_az, hstat_elev





