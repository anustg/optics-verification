import numpy as np
import pylab
np.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.sphere_surface import *
from tracer.paraboloid import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.models.one_sided_mirror import *
import types
import time

import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

'''
_________________________________________________________________________________________________________________
John.py:

Test on cavity geometry using 2 frustii.
_________________________________________________________________________________________________________________
'''
t0 = time.clock()

A = Assembly()

# Paraboloidal dish...
alpha = 0.1 # dish absorptivity
d = 22 # dish diameter
f = 13.4 # dish focal length
tr0 = translate(z=-f)
#P = AssembledObject(surfs=[Surface(ParabolicDishGM(d, f), Reflective(absorptivity=alpha))], transform=tr0)
P = AssembledObject(surfs=[Surface(ParabolicDishGM(d, f), RealReflective(absorptivity=alpha, sigma_xy=4e-3))], transform=tr0)
A.add_object(P)

# A beautiful thick, double frustum with 2 different sides. Still leaks.
alpha2 = 0.5
tr = N.dot(rotx(N.pi), translate(z=-0.3))
width = 0.001 # receiver thickness at frustii junction
# 1st frustum
CO1front = Surface(ConicalFrustum(z1=-0.5,r1=0.01,z2=0,r2=0.7), Reflective(alpha2))
CO1back = Surface(ConicalFrustum(z1=-0.5,r1=0.01+width,z2=0,r2=0.7), Reflective(alpha2))
CO1 = AssembledObject(surfs=[CO1front,CO1back], transform=tr)
CO1.surfaces_for_next_iteration = types.MethodType(surfaces_for_next_iteration, CO1, CO1.__class__)
# 2nd frustum
CO2front = Surface(ConicalFrustum(z1=0,r1=0.7,z2=0.3,r2=0.4), Reflective(alpha2))
CO2back = Surface(ConicalFrustum(z1=0,r1=0.7+width,z2=0.3,r2=0.4), Reflective(alpha2))
CO2 = AssembledObject(surfs=[CO2front,CO2back], transform=tr)
CO2.surfaces_for_next_iteration = types.MethodType(surfaces_for_next_iteration, CO2, CO2.__class__)
# add to general assembly
A.add_object(CO1)
A.add_object(CO2)

# A cylinder!
CY1 = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=1, height=1), AbsorberReflector(0.1))], transform=rotx(N.pi/2))
#A.add_object(CY1)
 
# Source definition
cr = np.array([[0,0,2*f]]).T
dr = np.array([0,0,-1])
ar = 0 # radians, sun rays angular range (what's the correct value?)
G = 1000. # W/m2 solar flux
nrays = 1000 # number of rays escaping the source
#TODO code in the Buie sunshape instead of a pillbox
src = solar_disk_bundle(nrays, cr, dr, d*0.6, ar, G)

# Raytrace!
engine = TracerEngine(A)
itmax = 1000 # stop iteration after this many ray bundles were generated (i.e. 
            # after the original rays intersected some surface this many times).
minener = 0.001 # minimum energy threshold
engine.ray_tracer(src, itmax, minener)

t1 = time.clock()-t0
print 'Raytrace calculation time: ', t1

'''
__________________________________________________________________________________________________________________
Rendering:

Renders the scene. Offers the option to highlight specific rays according to the number of times they have been 
reflected.
__________________________________________________________________________________________________________________
'''

def show_rays(engine, escaping_len=5., highlight_level=None):
    """
    Function to draw the rays to a Coin3D scenegraph.
    """
    tree = engine.tree
    no = coin.SoSeparator()
    #print tree.num_bunds()
    
    # loop through the reflection sequences?
    co = [] # regular lines
    co_h = [] # highlighted lines
    pos = [] # 2D level text position
    text = [] # 2D level text
    hist = {} # ray histories, for highlighted rays

    for level in xrange(tree.num_bunds()):
        start_rays = tree[level]
        sv = start_rays.get_vertices()
        sd = start_rays.get_directions()
        se = start_rays.get_energy()
        
        if level == tree.num_bunds() - 1:
            parents = []
        else:
            end_rays = tree[level + 1]
            ev = end_rays.get_vertices()
            parents = end_rays.get_parents()

        # loop through individual rays in this bundle
        for ray in xrange(start_rays.get_num_rays()):
            if se[ray] <= minener:
                # ignore rays with starting energy smaller than energy cutoff
                continue
            
            if ray in parents:
                # Has a hit on another surface
                first_child = N.where(ray == parents)[0][0]
                c1 = sv[:,ray]
                c2 = ev[:,first_child]
                #endpoints = N.c_[sv[:,ray], ev[:,first_child]]
            else:
                l = escaping_len
                if level == 0:
                    l = 0.1
                # Escaping ray.
                c1 = sv[:,ray]
                c2 = sv[:,ray] + sd[:,ray]*l
            if level == highlight_level:
                # Highlight rays that have the highlight level set.
                hist[ray] = tree.ray_history(ray,level)
                co_h += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]
            else:
                co += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]
            # Position and text of the level 2D text. ratio is the parameter of the text position on the ray.
            # eg. 1/5 is equivalent to a text position at one fifth of the total ray length in the scene.
            ratio = 1./5
            c3 = ratio*(((1./ratio)-1.)*c1+c2)
            pos += [(c3[0],c3[1],c3[2])]        
            text.append(str(level))

    def plot_rays_color(co, color=(1,1,0.5)):
        """
        Add ray set of line color `color` to scenegraph `node`. Rays `co`
        should be stored as sequences of 3-vector pairs in a list, eg
        [(x1,y1,z1),(x2,y2,z2),...]
        """
        no1 = coin.SoSeparator()

        ma1 = coin.SoMaterial()
        ma1.diffuseColor = color
        no1.addChild(ma1)

        ds = coin.SoDrawStyle()
        ds.style = ds.LINES
        ds.lineWidth = 2
        no1.addChild(ds)

        coor = coin.SoCoordinate3()
        coor.point.setValues(0, len(co), co)
        no1.addChild(coor)

        ls = coin.SoLineSet()
        ind = [2] * (len(co)/2)
        ls.numVertices.setValues(0, len(ind), ind)
        no1.addChild(ls)

        return no1

    def plot_level_number(text_pos, level_number):
        """
        Shows the number of reflections a ray has had as a 2D text on the scene.
        Arguments:
        text_pos - the position of the 2D text over the ray
        level_number - the number of reflections already encountered by the ray according to ray history.
        """
        no2 = coin.SoSeparator()
             
        tr = coin.SoTransform()
        tr.translation.setValue(text_pos)
        no2.addChild(tr)

        fo = coin.SoFont()
        fo.name.setValue("Arial-Bold")
        fo.size.setValue(15)
        no2.addChild(fo)   

        ma2 = coin.SoMaterial()
        ma2.diffuseColor.setValue(1,0,1)
        no2.addChild(ma2) 

        tx = coin.SoText2()      
        tx.string = level_number        
        no2.addChild(tx)
                
        return no2       
    
    no.addChild(plot_rays_color(co))
    no.addChild(plot_rays_color(co_h, color=(1,1,1)))
    num = len(text)    
    #for ref in range(num):
    #    no.addChild(plot_level_number(pos[ref], text[ref]))
    print "Number of reflections", num
    return no

def detailed_ray_history(engine,seq):
    """
    Print a detailed ray history for a particular ray, being a tuple returned 
    from TraceTree.ray_history.
    """
    tree = engine.tree
    n = len(seq)
    for i in range(n):
        bund = tree[i]
        ray = seq[(n-1)-i]
        print "...bund",i,"ray",ray
        print "...from",bund.get_vertices()[:,ray]
        print "...direction",bund.get_directions()[:,ray]

def axis_labels(length=1):
    """
    Create axis arrows/labels for addition to a Coin3D scenegraph.
    """
    r = coin.SoSeparator()
    st = coin.SoDrawStyle()
    r.addChild(st)
    st.lineWidth=3
    data = {'x':(1,0,0), 'y':(0,1,0), 'z':(0,0,1)}
    for k in data:
        vx,vy,vz = data[k]
        vec = (length*vx, length*vy, length*vz)
        
        s1 = coin.SoSeparator()

        la = coin.SoLabel()
        la.label = k
        s1.addChild(la)

        tr1 = coin.SoTranslation()
        tr1.translation = vec
        s1.addChild(tr1)

        r.addChild(s1)
    
        s2 = coin.SoSeparator()

        tr2 = coin.SoTransform()
        tr2.translation.setValue(data[k])
        s2.addChild(tr2)

        matxt = coin.SoMaterial()
        matxt.diffuseColor = data[k]
        s2.addChild(matxt)
       
        txaxis = coin.SoText2()      
        txaxis.string = k       
        s2.addChild(txaxis)

        r.addChild(s2)

        ma = coin.SoMaterial()
        ma.diffuseColor = data[k]
        r.addChild(ma)

        co = coin.SoCoordinate3()
        co.point.setValues(0,2,[(0,0,0),vec])
        r.addChild(co)

        ls = coin.SoLineSet()
        ls.numVertices.setValues(0,1,[2])
        r.addChild(ls)
    return r

# Render the scene with Pivy
r = coin.SoSeparator()

r.addChild(axis_labels())
r.addChild(show_rays(engine, highlight_level=None))
r.addChild(A.get_scene_graph())

win = SoGui.init("hello")
viewer = SoGuiExaminerViewer(win)
viewer.setSceneGraph(r)
viewer.setTitle("Examiner Viewer")
viewer.viewAll()
viewer.show()

SoGui.show(win)
SoGui.mainLoop()

# vim: et:ts=4
