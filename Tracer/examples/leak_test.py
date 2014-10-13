import numpy as np
import pylab
np.set_printoptions(linewidth=140)

from tracer.surface import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.trace_tree import *

import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

'''
_________________________________________________________________________________________________________________
leak_test:

Test of cylindrical geometry hit by a simple sources to figure out if some specific ray incidence causes the 
leak. All single_ray_sources are added to the srcs list and concatenated. 
_________________________________________________________________________________________________________________
'''

A = Assembly()
'''
<<< Geometry declaration >>>
'''
# Cylinder.
alpha1 = 0.8
CYL = AssembledObject(surfs=[Surface(FiniteCylinder(diameter=1,height=1), ReflectiveReceiver(alpha1))], transform = N.dot(rotx(N.pi/2.), translate(z=-0.5)))
A.add_object(CYL)
# Cone
CON = AssembledObject(surfs=[Surface(FiniteCone(r=1.,h=1.), ReflectiveReceiver(alpha1))], transform = N.dot(rotx(N.pi/2.), translate(x=5, z=-1)))
A.add_object(CON)
# Frustum
FRU = AssembledObject(surfs=[Surface(ConicalFrustum(z1=0,r1=1,z2=1,r2=0.5), ReflectiveReceiver(alpha1))], transform = translate(x=10.))
A.add_object(FRU)
'''
<<< Sources declaration >>>
'''
# Sources: define as many single ray sources as needed.
srcs = []
G = 1000. # W/m2 solar flux
'''
# Cylinder sources:
# Cylinder Source 1:
centerc1 = np.array([[0,3,-3]]).T
directc1 = np.array([0,-1.2,1])
srcc1 = single_ray_source(centerc1, directc1, G)
srcs.append(srcc1)
# Cylinder Source 2:
centerc2 = np.array([[-0.4,0.2,3]]).T
directc2 = np.array([0,0,-1])
srcc2 = single_ray_source(centerc2, directc2, G)
srcs.append(srcc2)
# Cylinder Source 3:
centerc3 = np.array([[3,2.5,0]]).T
directc3 = np.array([-1,-1,0])
srcc3 = single_ray_source(centerc3, directc3, G)
srcs.append(srcc3)
# Cylinder Source 4:
centerc4 = np.array([[3,2.7,0]]).T
directc4 = np.array([-1,-1,0])
srcc4 = single_ray_source(centerc4, directc4, G)
srcs.append(srcc4)
# Cylinder Source 5:
centerc5 = np.array([[3,2.9,0]]).T
directc5 = np.array([-1,-1,0])
srcc5 = single_ray_source(centerc5, directc5, G)
srcs.append(srcc5)
# Cylinder Source 6:
centerc6 = np.array([[3,3,0]]).T
directc6 = np.array([-1,-1,0])
srcc6 = single_ray_source(centerc6, directc6, G)
srcs.append(srcc6)
# Cylinder Source 7:
centerc7 = np.array([[-0.2,0.2,3]]).T
directc7 = np.array([0,0,-1])
srcc7 = single_ray_source(centerc7, directc7, G)
srcs.append(srcc7)
# Cylinder Source 8:
centerc8 = np.array([[0,0.2,3]]).T
directc8 = np.array([0,0,-1])
srcc8 = single_ray_source(centerc8, directc8, G)
srcs.append(srcc8)
'''  
# FiniteCone sources:
# FiniteCone Source 1:
centerco1 = np.array([[5,-3,0]]).T
directco1 = np.array([0,1,0])
srcco1 = single_ray_source(centerco1, directco1, G)
#srcs.append(srcco1)
# FiniteCone Source 2:
centerco2 = np.array([[5,0.3,3]]).T
directco2 = np.array([0,0,-1])
srcco2 = single_ray_source(centerco2, directco2, G)
srcs.append(srcco2)
# FiniteCone Source 3:
centerco3 = np.array([[5,0,0.2]]).T
directco3 = np.array([-1,0,0])
srcco3 = single_ray_source(centerco3, directco3, G)
srcs.append(srcco3)
# FiniteCone Source 4:
centerco4 = np.array([[5.25,3,0.2]]).T
directco4 = np.array([0,-1,0])
srcco4 = single_ray_source(centerco4, directco4, G)
srcs.append(srcco4) 
# FiniteCone Source 5:
centerco5 = np.array([[5,2,0]]).T
directco5 = np.array([0,-1,0])
srcco5 = single_ray_source(centerco5, directco5, G)
srcs.append(srcco5)
# FiniteCone Source 6:
centerco6 = np.array([[6,6,0]]).T
directco6 = np.array([-0.2,-1,0])
srcco6 = single_ray_source(centerco6, directco6, G)
srcs.append(srcco6)
# FiniteCone Source 7:
centerco7 = np.array([[4.6,0.3,3]]).T
directco7 = np.array([0,0,-1])
srcco7 = single_ray_source(centerco7, directco7, G)
srcs.append(srcco7)
# FiniteCone Source 8:
centerco8 = np.array([[5,-1,1]]).T
directco8 = np.array([0,1,-0.6])
srcco8 = single_ray_source(centerco8, directco8, G)
srcs.append(srcco8)
# FiniteCone Source 9:
centerco9 = np.array([[5,-1,1]]).T
directco9 = np.array([0,1,-0.7])
srcco9 = single_ray_source(centerco9, directco9, G)
srcs.append(srcco9)
# FiniteCone Source 10:
centerco10 = np.array([[5,-1,1]]).T
directco10 = np.array([0,1,-0.8])
srcco10 = single_ray_source(centerco10, directco10, G)
srcs.append(srcco10)
# FiniteCone Source 11:
centerco11 = np.array([[5,-1,1]]).T
directco11 = np.array([0,1,-0.9])
srcco11 = single_ray_source(centerco11, directco11, G)
srcs.append(srcco11)
 
# Frustum sources
# Frustum Source 1:
centerf1 = np.array([[7,0.2,0.2]]).T
directf1 = np.array([1,0,0])
srcf1 = single_ray_source(centerf1, directf1, G)
srcs.append(srcf1)
# Frustum Source 2:
centerf2 = np.array([[10,3,3]]).T
directf2 = np.array([0,-1,-1])
srcf2 = single_ray_source(centerf2, directf2, G)
srcs.append(srcf2)
# Frustum Source 3:
centerf3 = np.array([[10,0,0.2]]).T
directf3 = np.array([-1,0,0])
srcf3 = single_ray_source(centerf3, directf3, G)
srcs.append(srcf3)
# Frustum Source 4:
centerf4 = np.array([[10,0.6,3]]).T
directf4 = np.array([0,0,-1])
srcf4 = single_ray_source(centerf4, directf4, G)
srcs.append(srcf4)
# Frustum Source 5:
centerf5 = np.array([[7,0.8,0.2]]).T
directf5 = np.array([1,0,0])
srcf5 = single_ray_source(centerf5, directf5, G)
srcs.append(srcf5)  
# Frustum Source 6:
centerf6 = np.array([[7,0.6,0.2]]).T
directf6 = np.array([1,0,0])
srcf6 = single_ray_source(centerf6, directf6, G)
srcs.append(srcf6)  
# Frustum Source 7:
centerf7 = np.array([[7,0.4,0.2]]).T
directf7 = np.array([1,0,0])
srcf7 = single_ray_source(centerf7, directf7, G)
srcs.append(srcf7)  
# Frustum Source 8:
centerf8 = np.array([[10,0,-1]]).T
directf8 = np.array([0,0.5,1])
srcf8 = single_ray_source(centerf8, directf8, G)
srcs.append(srcf8)  
# Frustum Source 9:
centerf9 = np.array([[10,0,-1]]).T
directf9 = np.array([0,0.7,1])
srcf9 = single_ray_source(centerf9, directf9, G)
srcs.append(srcf9) 
# Frustum Source 10:
centerf10 = np.array([[10,0,-1]]).T
directf10 = np.array([0,0.8,1])
srcf10 = single_ray_source(centerf10, directf10, G)
srcs.append(srcf10)  
# Frustum Source 11:
centerf11 = np.array([[10,0,-1]]).T
directf11 = np.array([0,0.9,1])
srcf11 = single_ray_source(centerf11, directf11, G)
srcs.append(srcf11)   

# Concatenate Sources using ray_bundle method:
src = concatenate_rays(srcs)
'''
<<< RAYTRACE! >>>
'''
# TracerEngine instance declaration and ray_tracer parameter choice.
engine = TracerEngine(A)
itmax = 10 # stop iteration after this many ray bundles were generated (i.e. 
            # after the original rays intersected some surface this many times).
minener = 0.1 # minimum energy threshold
engine.ray_tracer(src, itmax, minener, tree = True)
Cone_acc = AbsorptionAccountant.get_all_hits(CON.get_surfaces()[0].get_optics_manager())
Cone_abs = Cone_acc[0]
Cone_hits = Cone_acc[1]
Frustum_acc = AbsorptionAccountant.get_all_hits(FRU.get_surfaces()[0].get_optics_manager())
Frustum_abs = Frustum_acc[0]
Frustum_hits = Frustum_acc[1]
print 'Cone_abs', Cone_abs
print 'Frustum_abs', Frustum_abs
print 'rays_in'
print engine.rays_in
print 'lost_energy'
print engine.lost_energy

print 'balance:', N.sum(src.get_energy())-N.sum(N.concatenate(engine.lost_energy))-N.sum(N.concatenate((Cone_abs,Frustum_abs)))
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
    for i in range(num):
        no.addChild(plot_level_number(pos[i], text[i]))
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
