import time
import numpy as N
from pathos.multiprocessing import ProcessingPool as Pool
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.assembly import *

class TracerEngineMP(TracerEngine):
	'''
	Famework for multi-processing using the tracer engine as is. 
	Requires pathos: https://github.com/uqfoundation/pathos

	Inheritance is broken by the multiprocessing pool and rebuilt on the tree and self._asm
	The original assembly needs to be reallocated after the simulation to be able to get the values stored in the optical managers with previously defined objects.

	Not the cleanest or finest implementation. Could be blended with the original engine and use the very same api. It works.
	'''

	def trace(self, source):
		# Raytrace in a separate method to be called by the ProcessingPool.map() method.
		self.ray_tracer(source, self.itmax, self.minener, tree=True)
		# Returns the tracer_engine instance traced.
		return self

	def multi_ray_sim(self, sources, procs=8):
		self.minener = 1e-10 # minimum energy threshold
		self.itmax = 1000 # stop iteration after this many ray bundles were generated (i.e. 
					# after the original rays intersected some surface this many times).
		# The multiprocessing raytracing method to call from the original engine.
		if len(sources) != procs:
			raise Exception('Number of sources and processors do not agree')

		# Creates a pool of processes and makes them raytrace one different source each. The resm list returned is a list of copies of the original engine post raytrace.
		pool = Pool(processes=procs)
		resm = pool.map(self.trace, sources)

		# New tree container and length envaluation to redimension it.
		tree_len = N.zeros(len(resm), dtype=N.int)
		trees = []

		for eng in xrange(len(resm)):
			# Get and regroup results in one tree and assembly only:
			S = resm[eng]._asm.get_surfaces()
			tree_len[eng] = len(resm[eng].tree._bunds)
			trees.append(resm[eng].tree)
			# Next loop is to get the optics callable objects and copy regroup their values without asumptions about what they are.
			for s in xrange(len(S)):
				part_res = S[s]._opt.__dict__
				keys = S[s]._opt.__dict__.keys()
				for k in xrange(len(keys)):
					if (keys[k] == '_opt') or (keys[k] == '_abs'):
						continue
					if len(self._asm.get_surfaces()[s]._opt.__dict__[keys[k]]) < 1:
						self._asm.get_surfaces()[s]._opt.__dict__[keys[k]] = part_res[keys[k]]
					elif len(part_res[keys[k]]) < 1:
						continue
					else:
						self._asm.get_surfaces()[s]._opt.__dict__[keys[k]][0] = N.append(self._asm.get_surfaces()[s]._opt.__dict__[keys[k]][0], part_res[keys[k]][0], axis=1)

		# Regroup trees:
		self.tree = RayTree() # Create a new tree for all
		for t in xrange(N.amax(tree_len)): # Browse through general tree levels up to the maximum length that has been raytraced
			for eng in xrange(len(resm)): # Browse through bundles of each parallel engine.
				if t<(tree_len[eng]): # to not go over the length of the present parallel tree.
					if t==len(self.tree._bunds): # if the index is greater than the actual length of the general tree, add a new bundle to the general tree with the present parallel bundle to initialise it.
						bundt = trees[eng]._bunds[t]
					else:	
						if t>0: # adapt parents indexing prior to concatenation
							trees[eng]._bunds[t].set_parents(trees[eng]._bunds[t].get_parents()+len(self.tree._bunds[t].get_parents()))
						bundt = concatenate_rays([bundt, trees[eng]._bunds[t]])
			self.tree.append(bundt)
		
		trees = 0
