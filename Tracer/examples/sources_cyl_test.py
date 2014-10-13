import numpy as N
import pylab
N.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.cone import *
from tracer.flat_surface import *
from tracer.cylinder import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.trace_tree import *
from tracer.CoIn_rendering.rendering import *

rc = 1.
lc = 1.

num_rays = 100000
el_CYL = 5
precision = 1e-3

VF = N.zeros((1+el_CYL+1, 1+el_CYL+1))
VF_avg = VF
test = VF
areas = N.zeros(N.shape(VF)[0])


def reset_opt():
	AP.get_surfaces()[0].get_optics_manager().reset()
	CYL.get_surfaces()[0].get_optics_manager().reset()
	BOT.get_surfaces()[0].get_optics_manager().reset()

# In final implementation use the object's assembly and add the aperture cover to perform VF calculations.
A = Assembly()

CYL = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=2.*rc,height=lc), ReflectiveReceiver(1.))], transform = translate(z=lc/2.))
AP = AssembledObject(surfs=[Surface(RoundPlateGM(Re=rc), ReflectiveReceiver(1.))], transform = None)
BOT = AssembledObject(surfs=[Surface(RoundPlateGM(Re=rc), ReflectiveReceiver(1.))], transform = N.dot(translate(z=lc),rotx(N.pi)))


A.add_object(CYL)
A.add_object(AP)
A.add_object(BOT)

engine = TracerEngine(A)

# Test stuff
engine.Scene = Assembly()
engine.Scene.add_assembly(A)

itmax = 1 # stop iteration after this many ray bundles were generated (i.e. 
            # after the original rays intersected some surface this many times).
engine.minener = 0.01/num_rays # minimum energy threshold

areas[0] = N.pi*rc**2

#for p in xrange(passes):
p=0
while 1:
	print 'PASS: ',p
	reset_opt()

	SA = solar_disk_bundle(num_rays, center=N.vstack([0,0,0]), direction=N.array([0,0,1]), radius=rc, ang_range=N.pi/2., flux=1./(N.pi*rc**2))
	engine.ray_tracer(SA, itmax, engine.minener, tree = True)

	Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
	Cylinder_abs, Cylinder_hits = CYL.get_surfaces()[0].get_optics_manager().get_all_hits()
	Bottom_abs, Bottom_hits = BOT.get_surfaces()[0].get_optics_manager().get_all_hits()
	
	for i in range(N.shape(VF)[0]):
		if i == 0:		
			VF[0,i] = N.sum(Aperture_abs)
		elif i < (el_CYL+1):
			VF[0,i] = N.sum(Cylinder_abs[N.where(N.logical_and(Cylinder_hits[2]>=((i-1)*lc/el_CYL),Cylinder_hits[2]<(i*lc/el_CYL)))])
		else:
			VF[0,i] = N.sum(Bottom_abs)

	for elc in range(el_CYL):

		reset_opt()
		areas[elc+1] = 2.*N.pi*rc*lc/el_CYL

		S = vf_cylinder_bundle(num_rays, rc=rc, lc=lc/el_CYL, center=N.vstack([0,0,elc*lc/el_CYL]), direction=N.array([0,0,1]), rays_in=True)
		engine.ray_tracer(S, itmax, engine.minener, tree = True)

		Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
		Cylinder_abs, Cylinder_hits = CYL.get_surfaces()[0].get_optics_manager().get_all_hits()
		Bottom_abs, Bottom_hits = BOT.get_surfaces()[0].get_optics_manager().get_all_hits()

		for i in range(N.shape(VF)[0]):
			if i == 0:		
				VF[1+elc,i] = N.sum(Aperture_abs)
			elif i < (el_CYL+1):
				VF[1+elc,i] = N.sum(Cylinder_abs[N.where(N.logical_and(Cylinder_hits[2]>=((i-1)*lc/el_CYL),Cylinder_hits[2]<(i*lc/el_CYL)))])
			else:
				VF[1+elc,i] = N.sum(Bottom_abs)

	reset_opt()
	areas[1+el_CYL] = N.pi*rc**2

	SB = solar_disk_bundle(num_rays, center=N.vstack([0,0,lc]), direction=N.array([0,0,-1]), radius=rc, ang_range=N.pi/2., flux=1./(N.pi*rc**2))
	engine.ray_tracer(SB, itmax, engine.minener, tree = True)

	Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
	Cylinder_abs, Cylinder_hits = CYL.get_surfaces()[0].get_optics_manager().get_all_hits()
	Bottom_abs, Bottom_hits = BOT.get_surfaces()[0].get_optics_manager().get_all_hits()

	for i in range(N.shape(VF)[0]):
		if i == 0:		
			VF[1+el_CYL,i] = N.sum(Aperture_abs)
		elif i < (el_CYL+1):
			VF[1+el_CYL,i] = N.sum(Cylinder_abs[N.where(N.logical_and(Cylinder_hits[2]>=((i-1)*lc/el_CYL),Cylinder_hits[2]<(i*lc/el_CYL)))])
		else:
			VF[1+el_CYL,i] = N.sum(Bottom_abs)

	VF_avg = (VF_avg*p+VF)/(p+1.)
	print 'VF'
	print VF
	print 'VF_avg'
	print VF_avg
	stdev_VF = N.sqrt(1./(p+1.)*(VF-VF_avg)**2)
	print 'stdev_VF'	
	print stdev_VF

	for i in range(N.shape(VF)[0]):
		test[i,:] = areas[i]*VF_avg[i,:]
	test = test-test.T

	print 'test'
	print test
	if p == 0:
		p+=1
	else:
		if N.logical_or((stdev_VF>precision).any(), (abs(test)/2.>precision).any()):
			p+=1
		else:
			break

print 'areas', areas
print 'VF_avg'
print N.sum(VF_avg, axis = 1)
print VF_avg





'''
case_render = Renderer(engine)
case_render.show_rays()


engine.ray_tracer(SA, itmax, engine.minener, tree = True)

case_render = Renderer(engine)
case_render.show_rays()
case.fluxmap()
plt.show()
'''
