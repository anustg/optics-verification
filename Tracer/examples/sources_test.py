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

ra = 0.5
rc = 1.
lf = 10.
lc = -9.1

num_rays = 1000
el_FRU = 1
el_CON = 1

def reset_opt():
	AP.get_surfaces()[0].get_optics_manager().reset()
	FRU.get_surfaces()[0].get_optics_manager().reset()
	CON.get_surfaces()[0].get_optics_manager().reset()

VF = N.zeros((1+el_FRU+el_CON, 1+el_FRU+el_CON))
test = N.zeros(N.shape(VF))
areas = N.zeros(1+el_FRU+el_CON)

passes = 1


# In final implementation use the object's assembly and add the aperture cover to perform VF calculations.
A = Assembly()

FRU = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0,r1=ra,z2=lf,r2=rc), ReflectiveReceiver(1.))], transform = None)
if lc<0.:
	CON = AssembledObject(surfs=[Surface(FiniteCone(r=rc, h=-lc), ReflectiveReceiver(1.))], transform=translate(z=lf+lc))
	rays_cone=False
else:
	CON = AssembledObject(surfs=[Surface(FiniteCone(r=rc, h=lc), ReflectiveReceiver(1.))], transform=N.dot(translate(z=lf+lc), rotx(N.pi)))
	rays_cone=True
AP = AssembledObject(surfs=[Surface(RoundPlateGM(Re=ra), ReflectiveReceiver(1.))], transform = None)

A.add_object(FRU)
A.add_object(CON)
A.add_object(AP)

engine = TracerEngine(A)

# Test stuff
engine.Scene = Assembly()
engine.Scene.add_assembly(A)

itmax = 1 # stop iteration after this many ray bundles were generated (i.e. 
            # after the original rays intersected some surface this many times).
engine.minener = 0.01/num_rays # minimum energy threshold

areas[0] = N.pi*ra**2

for p in xrange(passes):
	SA = solar_disk_bundle(num_rays, center=N.vstack([0,0,0]), direction=N.array([0,0,1]), radius=ra, ang_range=N.pi/2., flux=1./(N.pi*ra**2))
	engine.ray_tracer(SA, itmax, engine.minener, tree = True)

	Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
	Frustum_abs, Frustum_hits = FRU.get_surfaces()[0].get_optics_manager().get_all_hits()
	Cone_abs, Cone_hits = CON.get_surfaces()[0].get_optics_manager().get_all_hits()
	
for i in range(N.shape(VF)[0]):
	if i == 0:		
		VF[0,i] = N.sum(Aperture_abs)/passes
	elif i < (el_FRU+1):
		VF[0,i] = N.sum(Frustum_abs[N.where(N.logical_and(Frustum_hits[2]>=((i-1)*lf/el_FRU),Frustum_hits[2]<(i*lf/el_FRU)))])/passes
	else:
		VF[0,i] = N.sum(Cone_abs[N.where(N.logical_and(Cone_hits[2]>=(lf+(i-el_FRU-1)*lc/el_CON), Cone_hits[2]<(lf+(i-el_FRU)*lc/el_CON)))])/passes

for elf in range(el_FRU):

	reset_opt()

	areas[elf+1] = N.pi*(2.*ra+(2.*elf+1)*(rc-ra)/el_FRU)*N.sqrt((lf/el_FRU)**2+((rc-ra)/el_FRU)**2)

	for p in xrange(passes):

		S = vf_frustum_bundle(num_rays, r1=ra+elf*(rc-ra)/el_FRU, r2=ra+(elf+1)*(rc-ra)/el_FRU, depth=lf/el_FRU, center=N.vstack([0,0,elf*lf/el_FRU]), direction=N.array([0,0,1]), rays_in=True)
		engine.ray_tracer(S, itmax, engine.minener, tree = True)

		Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
		Frustum_abs, Frustum_hits = FRU.get_surfaces()[0].get_optics_manager().get_all_hits()
		Cone_abs, Cone_hits = CON.get_surfaces()[0].get_optics_manager().get_all_hits()

	for i in range(N.shape(VF)[0]):
		if i == 0:
			VF[1+elf,i] = N.sum(Aperture_abs)/passes
		elif i < (el_FRU+1):
			VF[1+elf,i] = N.sum(Frustum_abs[N.where(N.logical_and(Frustum_hits[2]>=((i-1)*lf/el_FRU),Frustum_hits[2]<(i*lf/el_FRU)))])/passes
		else:
			VF[1+elf,i] = N.sum(Cone_abs[N.where(N.logical_and(Cone_hits[2]>=(lf+(i-el_FRU-1)*lc/el_CON), Cone_hits[2]<(lf+(i-el_FRU)*lc/el_CON)))])/passes



for elc in range(el_CON):

	reset_opt()

	areas[elc+1+el_FRU] = N.pi*(2.*rc+(2.*elc+1)*(-rc)/el_CON)*N.sqrt((lc/el_CON)**2+((rc/el_CON)**2))

	for p in xrange(passes):

		S = vf_frustum_bundle(num_rays, r1=rc+elc*(-rc)/el_CON, r2=rc+(elc+1)*(-rc)/el_CON, depth=lc/el_CON, center=N.vstack([0,0,lf+lc*elc/el_CON]), direction=N.array([0,0,1]), rays_in=rays_cone)

		engine.ray_tracer(S, itmax, engine.minener, tree = True)

		Aperture_abs, Aperture_hits = AP.get_surfaces()[0].get_optics_manager().get_all_hits()
		Frustum_abs, Frustum_hits = FRU.get_surfaces()[0].get_optics_manager().get_all_hits()
		Cone_abs, Cone_hits = CON.get_surfaces()[0].get_optics_manager().get_all_hits()

	for i in range(N.shape(VF)[0]):
		if i == 0:
			VF[1+el_FRU+elc,i] = N.sum(Aperture_abs)/passes
		elif i < (el_FRU+1):
			VF[1+el_FRU+elc,i] = N.sum(Frustum_abs[N.where(N.logical_and(Frustum_hits[2]>=((i-1)*lf/el_FRU),Frustum_hits[2]<(i*lf/el_FRU)))])/passes
		else:
			VF[1+el_FRU+elc,i] = N.sum(Cone_abs[N.where(N.logical_and(Cone_hits[2]>=(lf+(i-el_FRU-1)*lc/el_CON), Cone_hits[2]<(lf+(i-el_FRU)*lc/el_CON)))])/passes
case_render = Renderer(engine)
case_render.show_rays()



print 'areas', areas
print 'VF'
print N.sum(VF, axis = 1)
print VF

for i in range(N.shape(VF)[0]):
	test[i,:] = areas[i]*VF[i,:]
print test
test = test-test.T

print 'test'
print test



'''
case_render = Renderer(engine)
case_render.show_rays()


engine.ray_tracer(SA, itmax, engine.minener, tree = True)

case_render = Renderer(engine)
case_render.show_rays()
case.fluxmap()
plt.show()
'''
