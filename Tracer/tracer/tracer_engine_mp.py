import time
import numpy as N
from pathos.multiprocessing import ProcessingPool as Pool
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.assembly import *
from copy import copy

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
		self.ray_tracer(source, self.reps, self.minener, tree=True)
		# Returns the tracer_engine instance traced.
		return self

	def multi_ray_sim(self, sources, procs=8):
		self.minener = 1e-10 # minimum energy threshold
		self.reps = 1000 # stop iteration after this many ray bundles were generated (i.e. 
					# after the original rays intersected some surface this many times).
		# The multiprocessing raytracing method to call from the original engine.
		if len(sources) != procs:
			raise Exception('Number of sources and processors do not agree')

		# Creates a pool of processes and makes them raytrace one different source each. The resm list returned is a list of copies of the original engine post raytrace.
		timetrace = time.clock()
		pool = Pool(processes=procs)
		resm = pool.map(self.trace, sources)

		timetrace = time.clock() - timetrace
		print 'Raytrace time: ', timetrace,'s'

		timepost = time.clock()

		# New general tree:
		for eng in xrange(len(resm)):
			# Get and regroup results in one tree and assembly only:
			if eng == 0 :
				self.tree =  resm[eng].tree
				self._asm = resm[eng]._asm
			else:
				eng_bunds = resm[eng].tree._bunds
				general_bunds = copy(self.tree._bunds)
				for b in xrange(len(eng_bunds)):
					if b>0:
						eng_bunds[b]._parents += general_bunds[b-1].get_num_rays()
					if b>=len(general_bunds):
						self.tree._bunds.append(eng_bunds[b])
					else:
						self.tree._bunds[b] = concatenate_rays([general_bunds[b], eng_bunds[b]])
				# Next loop is to get the optics callable objects and copy regroup their values without asumptions about what they are.
				surfs_engine = resm[eng]._asm.get_surfaces()
				surfs_general = self._asm.get_surfaces()

				for s in xrange(len(surfs_engine)):
					keys = surfs_engine[s]._opt.__dict__.keys()
					for k in keys:
						surf_engine_dict_k = surfs_engine[s]._opt.__dict__[k]
						surf_general_dict_k = surfs_general[s]._opt.__dict__[k]
						if k =='_opt' or k =='_abs' or k =='_sig':
							continue
						if len(surf_engine_dict_k):
							a1 = N.array(surf_engine_dict_k[0])
							if len(surf_general_dict_k):
								a2 = N.array(surf_general_dict_k[0])
								if a1.ndim < a2.ndim:
									a1.resize(a2.shape[0],1)
								if a1.ndim > a2.ndim:
									a2.resize(a1.shape[0],1)
								self._asm.get_surfaces()[s]._opt.__dict__[k] = [N.concatenate((a1,a2), axis=-1)]
							else:
								self._asm.get_surfaces()[s]._opt.__dict__[k] = [a1]

		timepost2 = time.clock()-timepost
		print 'Post processing reassociation time: ', timepost2,'s'
