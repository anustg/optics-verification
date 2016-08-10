import pivy.coin as coin
SOGUI_BINDING="SoQt"
from pivy.sogui import *

import numpy as N

class Renderer():
	'''	__________________________________________________________________________________________________
	Rendering:

	Renders the scene. Offers the option to highlight specific rays according to the number of times they have been 	reflected.	__________________________________________________________________________________________________

	'''
	def __init__(self, sim):
		self.sim = sim

		# Scene axis label
		length = 1
		self.r = coin.SoSeparator()
		st = coin.SoDrawStyle()
		st.lineWidth=3
		self.r.addChild(st)
		data = {'x':(1,0,0), 'y':(0,1,0), 'z':(0,0,1)}
		#bg = coin.SoGradientBackground()
		#bg.color0 = (1,1,1)
		#bg.color1 = (0,1,1)
		#self.r.addChild(bg)
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
		    self.r.addChild(s1)		

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
		    self.r.addChild(s2)

		    ma = coin.SoMaterial()
		    ma.diffuseColor = data[k]
		    self.r.addChild(ma)

		    co = coin.SoCoordinate3()
		    co.point.setValues(0,2,[(0,0,0),vec])
		    self.r.addChild(co)

		    ls = coin.SoLineSet()
		    ls.numVertices.setValues(0,1,[2])
		    self.r.addChild(ls)


	def show_geom(self, resolution=None):
		"""
		Method to draw the geometry of the scene to a Coin3D scenegraph.
		"""
		self.r.addChild(self.sim._asm.get_scene_graph(resolution))
		win = SoGui.init()
		viewer = SoGuiExaminerViewer(win)
		viewer.setSceneGraph(self.r)

		viewer.setTitle("Examiner Viewer")
		viewer.show()
		SoGui.show(win)

		SoGui.mainLoop()

	def show_rays(self, escaping_len=.2, max_rays=None, resolution=None):
		"""
		Method to draw the rays to a Coin3D scenegraph. Needs to be called after a raytrace has been peroformed.
		"""

		tree = self.sim.tree
		no = coin.SoSeparator()
		
		# loop through the reflection sequences?
		co = [] # regular lines
		pos = [] # 2D level text position
		text = [] # 2D level text
		hist = {} # ray histories, for highlighted rays

		for level in xrange(tree.num_bunds()):

			start_rays = tree[level]
			if max_rays ==  None:
				sv = start_rays.get_vertices()
				sd = start_rays.get_directions()
				se = start_rays.get_energy()
				rays = start_rays.get_num_rays()
			else:
				sv = start_rays.get_vertices(selector = N.arange(max_rays))
				sd = start_rays.get_directions(selector = N.arange(max_rays))
				se = start_rays.get_energy(selector = N.arange(max_rays))
				rays = max_rays

			if level == tree.num_bunds() - 1:
				parents = []
			else:
				end_rays = tree[level + 1]
				ev = end_rays.get_vertices()
				parents = end_rays.get_parents()

				if max_rays != None:
					in_max_rays = parents[:] <= max_rays					
					parents = N.argwhere(in_max_rays)
	
		    # loop through individual rays in this bundle
			for ray in xrange(rays):
				if se[ray] <= self.sim.minener:
					# ignore rays with starting energy smaller than energy cutoff
					continue
		
				if ray in parents:
					# Has a hit on another surface
					first_child = N.where(ray == parents)
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
				co += [(c1[0],c1[1],c1[2]), (c2[0],c2[1],c2[2])]

		color=(1,1,0.5)

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

		no.addChild(no1)
		self.r.addChild(no)
		self.r.addChild(self.sim._asm.get_scene_graph(resolution))

		win = SoGui.init()
		viewer = SoGuiExaminerViewer(win)
		viewer.setSceneGraph(self.r)
		viewer.setTitle("Examiner Viewer")

		viewer.show()
		SoGui.show(win)
		SoGui.mainLoop()
