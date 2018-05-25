import numpy as N
from numpy import random, linalg as LA
from tracer.ray_bundle import RayBundle
from tracer.spatial_geometry import *


def pillbox_sunshape_directions(num_rays, ang_rang):

    """
    Calculates directions for a ray bundles with ``num_rays`` rays, distributed
    as a Gaussian sunshape shining toward the +Z axis, and deviating from it by
    at most sigma (the standard deviation), such that if all rays have the same energy, the flux distribution comes out right.
    
    Arguments:
    num_rays - number of rays to generate directions for.
    ang_rang- in radians, the angular width of the solar disk
    
    Returns:
    A (3, num_rays) array whose each column is a unit direction vector for one
        ray, distributed to match a Gaussian sunshape.
    """

    x1=N.random.uniform(size=num_rays)
    x2=N.random.uniform(size=num_rays)

    theta=N.arcsin(N.sqrt(x1*(sin(ang_rang))**2))
    phi=x2*2.*N.pi

    a = N.vstack((N.cos(phi)*N.sin(theta), N.sin(phi)*N.sin(theta) , N.cos(theta)))

    return a

def pillbox_rect_bundle(num_rays, center, direction, x, y, ang_rang, flux=None):
    '''
    generate a rectancular ray bundle for pillbox sunshape

    Arguments:
    num_rays - number of rays to generate
    center - a column 3-array with the 3D coordinate of the ray bundle's center 
    direction - a 1D 3-array with the unit average direction vector for the bundle.
    x - width of the rectangular ray bundle
    y - height of the rectangular ray bundle
    ang_rang - the angular width of the solar disk
    flux - if not None, the ray bundle's energy is set such that each ray has
        an equal amount of energy, and the total energy is flux*pi*radius**2

    Returns:
    A RayBundle object with the above characteristics set.
    '''

    a= pillbox_sunshape_directions(num_rays, ang_rang)
    x*=1.2
    y*=1.2

    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)

	# Locations:
	# See [1]
    xs = N.random.uniform(low=-x/2., high=x/2., size=num_rays)
    ys = N.random.uniform(low=-y/2., high=y/2., size=num_rays)

    #if (direction == N.array([0,0,-1])).all():
    #    xs, ys = ys, xs

    # Rotate locations to the plane defined by <direction>:
    vertices_local = N.vstack((xs, ys, N.zeros(num_rays)))
    vertices_global = N.dot(perp_rot, vertices_local)

    #dirct=N.vstack((-N.ones(num_rays)*direction[0],-N.ones(num_rays)*direction[1],-N.ones(num_rays)*direction[2]))
    #vertices_global=vertices_local+dirct*100.


    rayb = RayBundle(vertices=vertices_global + center, directions=directions)



    if flux != None:
        rayb.set_energy(x*y/num_rays*flux*N.ones(num_rays))

    return rayb


def pillbox_effective_rays(num_rays, X,Y,Z,hits, direction, A_source, ang_rang,DNI):
    '''
    Generate a ray bundle according to Buie et al.: "Sunshape distributions for terrestrial simulations." Solar Energy 74 (2003) 113-122 (DOI: 10.1016/S0038-092X(03)00125-7).

    *the ray source just cover each individule mirror	

    Arguments:
    num_rays - number of rays over one heliostat, therefore total num of rays is num_rays*num_helios
    direction - (1,3)direction of the normal to the source disc.
    vertices - (3,n) array vertices of each ray
    energy - float - energy of each ray , W
    DNI - Direct normal irradiantion W/m2

    Returns:
    A raybundle object with the above characteristics set.
    '''

    total_rays=num_rays*N.sum(hits)

    a= pillbox_sunshape_directions(total_rays, ang_rang)
    # Rotate to a frame in which <direction> is Z:
    perp_rot = rotation_to_z(direction)	
    directions = N.sum(perp_rot[...,None] * a[None,...], axis=1)


    Xs=N.array([])
    Ys=N.array([])
    Zs=N.array([])

    #X,Y=N.meshgrid(X,Y)

    for i in xrange(len(hits)):
        for j in xrange(len(hits[i])):
            if hits[i,j]==1:
                xs = N.random.uniform(low=X[i], high=X[i+1], size=num_rays)
                ys = N.random.uniform(low=Y[j], high=Y[j+1], size=num_rays)  
                zs= Z*N.ones(num_rays)
                Xs=N.append(Xs,xs)  
                Ys=N.append(Ys,ys) 
                Zs=N.append(Zs,zs) 


    vertices=N.vstack((Xs,Ys,Zs))            


    energy=N.ones(total_rays)*DNI*A_source/float(total_rays)

    rayb = RayBundle(vertices = vertices, directions = directions)
    if DNI != None:
        rayb.set_energy(energy)

    print total_rays


    return rayb
       


