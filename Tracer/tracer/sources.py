"""
This module contains functions to create some frequently-used light sources.
Each function returns a RayBundle instance that represents the distribution of
rays expected from that source.

References:
.. [1] Monte Carlo Ray Tracing, Siggraph 2003 Course 44


Ye changed def(solar_rect_bundel)
Line189 vertices_local = N.vstack((ys, xs, N.zeros(num_rays)))
to      vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
28 Oct 2016
"""

from numpy import random, linalg as LA
import numpy as N
from tracer.ray_bundle import RayBundle
from tracer.spatial_geometry import *

def single_ray_source(position, direction, flux=None):
    '''
    Establishes a single ray source originating from a definned point on a defined exact 
    direction for the purpose of testing single ray behviours.

    Arguments:
    position - column 3-array with the ray's starting position.
    direction - a 1D 3-array with the unit average direction vector for the
                bundle.
    flux - if not None, the energy transported by the ray.

    Returns:
    A Raybundle object with the corresponding characteristics.
    '''
    directions = N.tile(direction[:,None],1)
    directions /= N.sqrt(N.sum(directions**2, axis=0))
    singray = RayBundle(vertices = position, directions = directions)
    singray.set_energy(flux*N.ones(1))
    return singray

def pillbox_sunshape_directions(num_rays, ang_range):
    """
    Calculates directions for a ray bundles with ``num_rays`` rays, distributed
    as a pillbox sunshape shining toward the +Z axis, and deviating from it by
    at most ang_range, such that if all rays have the same energy, the flux
    distribution comes out right.
    
    Arguments:
    num_rays - number of rays to generate directions for.
    ang_range - in radians, the maximum deviation from +Z.
    
    Returns:
    A (3, num_rays) array whose each column is a unit direction vector for one
        ray, distributed to match a pillbox sunshape.
    """
    # Diffuse divergence from +Z:
    # development based on eq. 2.12  from [1]
    xi1 = random.uniform(high=2.*N.pi, size=num_rays)
    xi2 = random.uniform(size=num_rays)
    theta = N.arcsin(N.sin(ang_range)*N.sqrt(xi2))

    sin_th = N.sin(theta)
    a = N.vstack((N.cos(xi1)*sin_th, N.sin(xi1)*sin_th , N.cos(theta)))

    return a

def bivariate_directions(num_rays, ang_range_hor, ang_range_vert):
    """
    Calculates directions for a ray bundles with ``num_rays`` rays, distributed
    as uniform bi-variate distribution shining toward the +Z axis and deviating from it by
    at most ang_range_hor on the zx plane and ang_range_vert on the yz plane such that if all rays have the same energy, the flux
    distribution comes out right.
    
    Arguments:
    num_rays - number of rays to generate directions for.
    ang_range_hor - in radians, the maximum deviation from +Z on the zx plane.
    ang_range_vert - in radians, the maximum deviation from +Z on the yz plane.
    
    Returns:
    A (3, num_rays) array whose each column is a unit direction vector for one
        ray, distributed to match a uniform bi-variate distribution.
    """
    # Diffuse divergence from +Z:
    # development based on eq. 2.12  from [1]
    '''
    xi1 = N.random.uniform(low=-1., high=1., size=num_rays)
    xi2 = N.random.uniform(low=-1., high=1., size=num_rays)
    
    theta_hor = N.arcsin(N.sin(ang_range_hor)*N.sqrt(xi1))
    theta_vert = N.arcsin(N.sin(ang_range_vert)*N.sqrt(xi2))

    xa = N.sin(theta_hor)
    ya = N.sin(theta_vert)
    za = N.sqrt(1.-N.sin(theta_hor)**2.-N.sin(theta_vert)**2.)

    a = N.vstack((xa, ya, za))
    '''
    return a

def edge_rays_directions(num_rays, ang_range):
    """
    Calculates directions for a ray bundles with ``num_rays`` rays, distributed
    as a pillbox sunshape shining toward the +Z axis, and deviating from it by
    at most ang_range, such that if all rays have the same energy, the flux
    distribution comes out right.
    
    Arguments:
    num_rays - number of rays to generate directions for.
    ang_range - in radians, the maximum deviation from +Z.
    
    Returns:
    A (3, num_rays) array whose each column is a unit direction vector for one
        ray, distributed to match a pillbox sunshape.
    """
    # Diffuse divergence from +Z:
    # development based on eq. 2.12  from [1]
    xi1 = random.uniform(high=2.*N.pi, size=num_rays)
    sin_th = N.ones(num_rays)*N.sin(ang_range)
    a = N.vstack((N.cos(xi1)*sin_th, N.sin(xi1)*sin_th , N.cos(N.ones(num_rays)*ang_range)))

    return a

def solar_disk_bundle(num_rays,  center,  direction,  radius, ang_range, flux=None, radius_in=0., angular_span=[0.,2.*N.pi], procs=1):
    """
    Generates a ray bundle emanating from a disk, with each surface element of 
    the disk having the same ray density. The rays all point at directions uniformly 
    distributed between a given angle range from a given direction.
    Setting of the bundle's energy is left to the caller.
    
    Arguments:
    num_rays - number of rays to generate.
    center - a column 3-array with the 3D coordinate of the disk's center
    direction - a 1D 3-array with the unit average direction vector for the
        bundle.
    radius - of the disk.
    ang_range - in radians, the maximum deviation from <direction>.
    flux - if not None, the ray bundle's energy is set such that each ray has
        an equal amount of energy, and the total energy is flux*pi*radius**2
    radius_in - Inner radius if the disc is pierced
    angular_span - wedge of the disc to consider
    
    Returns: 
    A RayBundle object with the above characteristics set.
    """

	# FIXME why should 'center' be a column vector... that's just annoying.

    radius = float(radius)
    radius_in = float(radius_in)
    a = pillbox_sunshape_directions(num_rays, ang_range)
        
    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)
    # Locations:
    # See [1]
    xi1 = random.uniform(size=num_rays)
    thetas = random.uniform(low=angular_span[0], high=angular_span[1], size=num_rays)
    rs = N.sqrt(radius_in**2.+xi1*(radius**2.-radius_in**2.))
    xs = rs * N.cos(thetas)
    ys = rs * N.sin(thetas)

    # Rotate locations to the plane defined by <direction>:
    vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
    vertices_global = N.dot(perp_rot, vertices_local)

    rayb = RayBundle(vertices=vertices_global + center, directions=directions)
    if flux != None:
        rayb.set_energy(N.pi*(radius**2.-radius_in**2.)/num_rays*flux*N.ones(num_rays))
    else:
        rayb.set_energy(N.ones(num_rays)/num_rays/procs)
    return rayb

def solar_rect_bundle(num_rays, center, direction, x, y, ang_range, flux=None):

    a = pillbox_sunshape_directions(num_rays, ang_range)

    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)

    xs = random.uniform(low=-x/2., high=x/2., size=num_rays)
    ys = random.uniform(low=-y/2., high=y/2., size=num_rays)

    if (direction == N.array([0,0,-1])).all():
        xs, ys = ys, xs

    # Rotate locations to the plane defined by <direction>:
    vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
    vertices_global = N.dot(perp_rot, vertices_local)

    rayb = RayBundle(vertices=vertices_global + center, directions=directions)
    if flux != None:
        rayb.set_energy(x*y/num_rays*flux*N.ones(num_rays))
    return rayb

#def bivariate_rect_bundle(num_rays, center, direction, x, y, ang_range_vert, ang_range_hor, flux=None):


def edge_rays_bundle(num_rays,  center,  direction,  radius, ang_range, flux=None, radius_in=0.):

    radius = float(radius)
    radius_in = float(radius_in)
    a = edge_rays_directions(num_rays, ang_range)
        
    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)
    # Locations:
    # See [1]
    xi1 = random.uniform(size=num_rays)
    thetas = random.uniform(high=2.*N.pi, size=num_rays)
    rs = N.sqrt(radius_in**2.+xi1*(radius**2.-radius_in**2.))
    xs = rs * N.cos(thetas)
    ys = rs * N.sin(thetas)

    # Rotate locations to the plane defined by <direction>:
    vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
    vertices_global = N.dot(perp_rot, vertices_local)

    rayb = RayBundle(vertices=vertices_global + center, directions=directions)
    if flux != None:
        rayb.set_energy(N.pi*(radius**2.-radius_in**2.)/num_rays*flux*N.ones(num_rays))
    return rayb

def buie_sunshape(num_rays, center, direction, radius, CSR, flux=None, pre_process_CSR=True):
    '''
    Generate a ray bundle according to Buie et al.: "Sunshape distributions for terrestrial simulations." Solar Energy 74 (2003) 113-122 (DOI: 10.1016/S0038-092X(03)00125-7).

    Arguments:
    num_rays - number of rays in the bundle
    center - position of the source center
    direction - direction of the normal to the source disc.
    radius - radius of the source disc
    CSR - Circumsolar ratio, fraction of the incoming solar energy which incident angle is greater than the angle subtended by the solar disc.
    flux - horizontal radiation density in W/m2

    Returns:
    A raybundle object with the above characteristics set.
    '''
    # Angles of importance:
    theta_dni = 4.65e-3 # mrad
    theta_tot = 43.6e-3 # mrad

    # Rays vertices (start positions):

    xv1 = random.uniform(size=num_rays)
    phiv = random.uniform(high=2.*N.pi, size=num_rays)
    rs = radius*N.sqrt(xv1)
    xs = rs * N.cos(phiv)
    ys = rs * N.sin(phiv)

    # Source surface area:
    S = N.pi*radius**2.

    # Uniform ray energy:
    energy = N.ones(num_rays)*flux*S/num_rays

    # Polar angle array:
    thetas = N.zeros(num_rays)

    # Discrete random ray directions generation according to Buie sunshape
    # Step 1: integration over the whole Sunshape: 
    nelem = 210

    theta_int = N.linspace(0., theta_dni*1000., num=nelem)
    phi_dni_int = N.sin(theta_int/1000.)*N.cos(0.326*theta_int)/N.cos(0.308*theta_int)
    integ_phi_dni = theta_dni/nelem/2.*(phi_dni_int[:-1]+phi_dni_int[1:])

    if CSR == 0.:
        integ_phi = N.sum(integ_phi_dni)
    else:
        if pre_process_CSR:
            if CSR<=0.1:
                CSR = -2.245e+03*CSR**4.+5.207e+02*CSR**3.-3.939e+01*CSR**2.+1.891e+00*CSR+8e-03
            else:
                CSR = 1.973*CSR**4.-2.481*CSR**3.+0.607*CSR**2.+1.151*CSR-0.020
        # Buie Sunshape parameters:
        kappa = 0.9*N.log(13.5*CSR)*CSR**(-0.3)
        gamma = 2.2*N.log(0.52*CSR)*CSR**(0.43)-0.1
        integ_phi_csr = 1e-6*N.exp(kappa)/(gamma+2.)*((theta_tot*1000.)**(gamma+2.)-(theta_dni*1000.)**(gamma+2.))
        integ_phi = N.sum(integ_phi_dni)+integ_phi_csr

    # Step 2: PDF and random variate declaration
    integ_pdf_dni = integ_phi_dni/integ_phi
    R_thetas = N.random.uniform(size=num_rays)

    # Step 3: polar angle determination: 
    aureole = R_thetas>=N.sum(integ_pdf_dni)
    for i in xrange(len(integ_pdf_dni)-1):
        dni_slice = N.logical_and((R_thetas >= N.sum(integ_pdf_dni[:i])), (R_thetas < N.sum(integ_pdf_dni[:i+1])))
        thetas[dni_slice] = theta_int[i]/1000.+N.random.uniform(size=N.sum(dni_slice))*theta_dni/400.

    if CSR>0.:
        thetas[aureole] = ((R_thetas[aureole]-1.)*((gamma+2.)/(10.**(3.*gamma)*N.exp(kappa))*N.sum(integ_phi_dni)-theta_dni**(gamma+2.))+R_thetas[aureole]*theta_tot**(gamma+2.))**(1./(gamma+2.))

    # Generate directions:
    xi1 = random.uniform(high=2.*N.pi, size=num_rays)
    sin_th = N.sin(N.hstack(thetas))
    a = N.vstack((N.cos(xi1)*sin_th, N.sin(xi1)*sin_th , N.cos(thetas)))

    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)

    # Rotate locations to the plane defined by <direction>:
    vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
    vertices_global = N.dot(perp_rot, vertices_local)
    
    rayb = RayBundle(vertices = vertices_global+center, directions = directions, energy = energy)

    return rayb


def square_bundle(num_rays, center, direction, width):
    """
    Generate a ray bundles whose rays are equally spaced along a square grid,
    and all pointing in the same direction.
    
    Arguments:
    num_rays - number of rays to generate.
    center - a column 3-array with the 3D coordinate of the disk's center
    direction - a 1D 3-array with the unit direction vector for the bundle.
    width - of the square of starting points.
    
    Returns: 
    A RayBundle object with the above charachteristics set.
    """
    rot = rotation_to_z(direction)
    directions = N.tile(direction[:,None], (1, num_rays))
    range1 = N.s_[-width:width:float(2*width)/N.sqrt(num_rays)]
    xs, ys = N.mgrid[range1, range1]
    vertices_local = N.array([xs.flatten(),  ys.flatten(),  N.zeros(len(xs.flatten()))])
    vertices_global = N.dot(rot,  vertices_local)

    rayb = RayBundle()
    rayb.set_vertices(vertices_global + center)
    rayb.set_directions(directions)
    return rayb

def vf_frustum_bundle(num_rays, r0, r1, depth, center, direction, flux=None , rays_in=True, procs=1, angular_span=[0.,2.*N.pi], angular_range=N.pi/2.):
    '''
    Generate a frustum shaped lambertian source with randomly situated rays to compute view factors. The overall energy of the bundle is 1.

    Arguments:
    num_rays - number of rays to generate.
    center - a column 3-array with the 3D coordinate of the center of one of the bases.
    r0 - The radius of the frustum base which center coordinate has been given.
    r1 - the radius of the frustum at the other base location.
    depth - the depth of the frustum.
    direction - The orientation of the overall bundle. 
               Positive if in the same direction as the depth.
    rays_in - True if rays are fired towards the axis of the frustum.
    angular_span - wedge of the shape to consider.

    Returns:
    A raybundle object with the above characteristics set.
    '''
    r0 = float(r0)
    r1 = float(r1)
    depth = float(depth)
    if r0>r1:
         raise AttributeError('wrong radii') 
    if depth <0.:
         raise AttributeError('wrong depth')
    num_rays = float(num_rays)

    dir_flat = pillbox_sunshape_directions(num_rays, angular_range)

    c = (r1-r0)/depth

    R = random.uniform(size=num_rays)

    zs = (-r0+N.sqrt(r0**2.+R*(r1**2.-r0**2.)))/c

    phi_s = random.uniform(low=angular_span[0], high=angular_span[1], size=num_rays)
    rs = r0+c*zs
    xs = rs * N.cos(phi_s)
    ys = rs * N.sin(phi_s)

    theta_s = N.arctan(c)
    theta_rot = -N.pi/2.+theta_s
    yrot = roty(theta_rot)[:3,:3]
    local_unit = N.zeros((N.shape(dir_flat)))
    for t in xrange(N.shape(dir_flat)[1]):
        zrot = rotz(phi_s[t])[:3,:3]
        rot = N.dot(zrot, yrot)
        local_unit[:,t] = N.dot(rot, dir_flat[:,t])

    if rays_in == False:
        local_unit = -local_unit

    vertices_local = N.vstack((xs, ys, zs))

    perp_rot = rotation_to_z(direction)
    vertices_global = N.dot(perp_rot, vertices_local)
    directions = N.dot(perp_rot, local_unit)

    if flux == None:
        energy = N.ones(num_rays)/num_rays/procs
    else:
        area = (angular_span[1]-angular_span[0])*(r1+r0)/2.*N.sqrt(abs(r1-r0)**2.+depth**2.)
        energy = N.ones(num_rays)*flux*area/num_rays/procs

    rayb = RayBundle(vertices = vertices_global+center, directions = directions, energy = energy)

    return rayb

def vf_cylinder_bundle(num_rays, rc, lc, center, direction, flux=None, rays_in=True, procs=1, angular_span=[0.,2.*N.pi]):
    '''
    Generate a cylinder shaped lambertian source with randomly situated rays to compute view factors. The overall energy of the bundle is 1.

    Arguments:
    num_rays - number of rays to generate.
    center - a column 3-array with the 3D coordinate of the center of one of the bases.
    rc - The radius of the cylinder.
    lc - the length of the cylinder.
    direction - the direction of outgoing rays as projected on the cylinder axis. 
               Positive if in the same direction as lc.
    rays_in - True if rays are fired towards the axis of the frustum.
    angular_span - wedge of the shape to consider.

    Returns:
    A raybundle object with the above characteristics set.
    '''
    rc = float(rc)
    lc = float(lc)
    num_rays = float(num_rays)

    zs = lc*random.uniform(size=num_rays)

    phi_s = random.uniform(low=angular_span[0], high=angular_span[1], size=num_rays)

    xs = rc * N.cos(phi_s)
    ys = rc * N.sin(phi_s)

    dir_flat = pillbox_sunshape_directions(num_rays, N.pi/2.)

    yrot = roty(-N.pi/2.)[:3,:3]
    local_unit = N.zeros((N.shape(dir_flat)))
    for t in range(N.shape(dir_flat)[1]):
        zrot = rotz(phi_s[t])[:3,:3]
        rot = N.dot(zrot, yrot)
        local_unit[:,t] = N.dot(rot, dir_flat[:,t])

    if rays_in == False:
        local_unit = -local_unit

    vertices_local = N.vstack((xs, ys, zs))
    perp_rot = rotation_to_z(direction)
    vertices_global = N.dot(perp_rot, vertices_local)
    directions = N.dot(perp_rot, local_unit)
    '''
    plt.hist(vertices_local[2,:]/(N.sqrt(vertices_local[0,:]**2.+vertices_local[1,:]**2.)))
    plt.show()
    '''
    if flux == None:
        energy = N.ones(num_rays)/num_rays/procs
    else:
        area = rc*(angular_span[1]-angular_span[0])*lc
        energy = N.ones(num_rays)*flux*area/num_rays/procs

    rayb = RayBundle(vertices = vertices_global+center, directions = directions, energy = energy)

    return rayb


# vim: et:ts=4
