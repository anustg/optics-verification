# -*- coding: utf-8 -*-
# A collection of callables and tools for creating them, that may be used for
# the optics-callable part of a Surface object.

from . import optics, ray_bundle, sources
from .spatial_geometry import rotation_to_z
import numpy as N

class Reflective(object):
    """
    Generates a function that represents the optics of an opaque, absorptive
    surface with specular reflections.
    
    Arguments:
    absorptivity - the amount of energy absorbed before reflection.
    
    Returns:
    refractive - a function with the signature required by Surface.
    """
    def __init__(self, absorptivity):
        self._abs = absorptivity
    
    def __call__(self, geometry, rays, selector):
        outg = rays.inherit(selector,
            vertices=geometry.get_intersection_points_global(),
            direction=optics.reflections(rays.get_directions()[:,selector], geometry.get_normals()),
            energy=rays.get_energy()[selector]*(1 - self._abs),
            parents=selector)
        return outg

class RealReflective(object):
    '''
    Generates a function that represents the optics of an opaque absorptive surface with specular reflections and realistic surface slope error. The surface slope error is considered equal in both x and y directions. The consequent distribution of standard deviation is described by a radial bivariate normal distribution law.

    Arguments:
    absorptivity - the amount of energy absorbed before reflection
    sigma_xy - Standard deviation of the reflected ray in the local x and y directions. 
    
    Returns:
    Reflective - a function with the signature required by surface
    '''
    def __init__(self, absorptivity, sigma_xy):
        self._abs = absorptivity
        self._sig = sigma_xy

    def __call__(self, geometry, rays, selector):
        ideal_normals = geometry.get_normals()
        # Creates projection of error_normal on the surface (sin can be avoided because of very small angles).
        normal_errors_x = N.sin(N.random.normal(scale=self._sig, size=N.shape(ideal_normals[1])))
        normal_errors_y = N.sin(N.random.normal(scale=self._sig, size=N.shape(ideal_normals[1])))
        normal_errors_z = N.zeros(N.shape(ideal_normals[1]))
        # Build the normal_error vectors in the local frame.
        normal_errors = N.dot(geometry._working_frame[:3,:3], N.vstack((normal_errors_x, normal_errors_y, normal_errors_z)))
        real_normals = ideal_normals + normal_errors
        real_normals_unit = real_normals/N.sqrt(N.sum(real_normals**2, axis=0))
        # Call reflective optics with the new set of normals to get reflections affected by 
        # shape error.
        outg = rays.inherit(selector,
            vertices = geometry.get_intersection_points_global(),
            direction = optics.reflections(rays.get_directions()[:,selector], real_normals_unit),
            energy = rays.get_energy()[selector]*(1 - self._abs),
            parents = selector)
        return outg

class AbsorptionAccountant(object):
    """
    This optics manager remembers all of the locations where rays hit it
    in all iterations, and the energy absorbed from each ray.
    """
    def __init__(self, real_optics, absorptivity, sigma_xy=None):
        """
        Arguments:
        real_optics - the optics manager class to actually use. Expected to
            have the _abs protected attribute, and accept absorptivity as its
            only constructor argument (as in Reflective and
            LambertianReflector below).
        absorptivity - to be passed to a new real_optics object.
        """
        if sigma_xy==None:
            self._opt = real_optics(absorptivity)
        else:
            self._opt = real_optics(absorptivity, sigma_xy)
        self.reset()
    
    def reset(self):
        """Clear the memory of hits (best done before a new trace)."""
        self._absorbed = []
        self._hits = []
    
    def __call__(self, geometry, rays, selector):
        self._absorbed.append(rays.get_energy()[selector]*self._opt._abs)
        self._hits.append(geometry.get_intersection_points_global())
        return self._opt(geometry, rays, selector)
    
    def get_all_hits(self):
        """
        Aggregate all hits from all stages of tracing into joined arrays.
        
        Returns:
        absorbed - the energy absorbed by each hit-point
        hits - the corresponding global coordinates for each hit-point.
        """
        if not len(self._absorbed):
            return N.array([]), N.array([]).reshape(3,0)
        
        return N.hstack([a for a in self._absorbed if len(a)]), \
            N.hstack([h for h in self._hits if h.shape[1]])

class ReflectiveReceiver(AbsorptionAccountant):
    """A wrapper around AbsorptionAccountant with a Reflective optics"""
    def __init__(self, absorptivity):
        AbsorptionAccountant.__init__(self, Reflective, absorptivity)

class RealReflectiveReceiver(AbsorptionAccountant):
    """A wrapper around AbsorptionAccountant with a RealReflective optics"""
    def __init__(self, absorptivity=0, sigma_xy=0):
        AbsorptionAccountant.__init__(self, RealReflective, absorptivity, sigma_xy)
        
class AbsorberReflector(Reflective):
    """
    This optics manager behaves similarly to the ReflectiveReceiver class,
    but adds directionality. In this way a simple one-side receiver doesn't
    necessitate an extra surface in the back.
    """
    def __call__(self, geometry, rays, selector):
        """
        Rays coming from the "up" side are reflected like in a Reflective
        instance, rays coming from the "down" side have their energy set to 0.
        As usual, "up" is the surface's Z axis.
        """
        outg = Reflective.__call__(self, geometry, rays, selector)
        energy = outg.get_energy()
        proj = N.sum(rays.get_directions()[:,selector] * geometry.up()[:,None], axis=0)
        energy[proj > 0] = 0
        outg.set_energy(energy)
        return outg

# added
class RealReflective_OneSide(RealReflective):
    """
    Adds directionality to an optics manager that is modelled to represent the
    optics of an opaque absorptive surface with specular reflections and realistic
    surface slope error.
    """
    def __call__(self, geometry, rays, selector):
        outg = RealReflective.__call__(self, geometry, rays, selector)
        energy = outg.set_energy()
        proj = N.sum(rays.get_directions()[:,selector]*geometry.up()[:,None], axis = 0)
        energy[proj > 0] = 0 # projection up - set energy to zero
        outg.set_energy(energy) #substitute previous step into ray energy array
        return outg
# added end
        
class RefractiveHomogenous(object):
    """
    Represents the optics of a surface bordering homogenous media with 
    constant refractive index on each side. The specific index in which a
    refracted ray moves is determined by toggling between the two possible
    indices.
    """
    def __init__(self, n1, n2):
        """
        Arguments:
        n1, n2 - scalars representing the homogenous refractive index on each
            side of the surface (order doesn't matter).
        """
        self._ref_idxs = (n1, n2)
    
    def toggle_ref_idx(self, current):
        """
        Determines which refractive index to use based on the refractive index
        rays are currently travelling through.

        Arguments:
        current - an array of the refractive indices of the materials each of 
            the rays in a ray bundle is travelling through.
        
        Returns:
        An array of length(n) with the index to use for each ray.
        """
        return N.where(current == self._ref_idxs[0], 
            self._ref_idxs[1], self._ref_idxs[0])
    
    def __call__(self, geometry, rays, selector):
        if len(selector) == 0:
            return ray_bundle.empty_bund()
        
        n1 = rays.get_ref_index()[selector]
        n2 = self.toggle_ref_idx(n1)
        refr, out_dirs = optics.refractions(n1, n2, \
            rays.get_directions()[:,selector], geometry.get_normals())
        
        if not refr.any():
            return perfect_mirror(geometry, rays, selector)
        
        # Reflected energy:
        R = N.ones(len(selector))
        R[refr] = optics.fresnel(rays.get_directions()[:,selector][:,refr],
            geometry.get_normals()[:,refr], n1[refr], n2[refr])
        
        # The output bundle is generated by stacking together the reflected and
        # refracted rays in that order.
        inters = geometry.get_intersection_points_global()
        reflected_rays = rays.inherit(selector, vertices=inters,
            direction=optics.reflections(
                rays.get_directions()[:,selector],
                geometry.get_normals()),
            energy=rays.get_energy()[selector]*R,
            parents=selector)
        
        refracted_rays = rays.inherit(selector[refr], vertices=inters[:,refr],
            direction=out_dirs, parents=selector[refr],
            energy=rays.get_energy()[selector][refr]*(1 - R[refr]),
            ref_index=n2[refr])
        
        return reflected_rays + refracted_rays

class LambertianReflector(object):
    """
    Represents the optics of an ideal diffuse (lambertian) surface, i.e. one
    that reflects rays in a random direction (uniform distribution of
    directions in 3D, see tracer.sources.pillbox_sunshape_directions)
    """
    def __init__(self, absorptivity):
        self._abs = absorptivity
    
    def __call__(self, geometry, rays, selector):
        """
        Arguments:
        geometry - a GeometryManager which knows about surface normals, hit
            points etc.
        rays - the incoming ray bundle (all of it, not just rays hitting this
            surface)
        selector - indices into ``rays`` of the hitting rays.
        """
        directs = sources.pillbox_sunshape_directions(len(selector), N.pi/2.)
        directs = N.sum(rotation_to_z(geometry.get_normals().T) * directs.T[:,None,:], axis=2).T
        
        outg = rays.inherit(selector,
            vertices=geometry.get_intersection_points_global(),
            energy=rays.get_energy()[selector]*(1 - self._abs),
            direction=directs, 
            parents=selector)
        return outg

class LambertianReceiver(AbsorptionAccountant):
    """A wrapper around AbsorptionAccountant with LambertianReflector optics"""
    def __init__(self, absorptivity):
        AbsorptionAccountant.__init__(self, LambertianReflector, absorptivity)


# vim: et:ts=4
