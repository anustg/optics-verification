import numpy as N
import matplotlib.pyplot as plt

from tracer.models.four_parameters_cavity import *
from tracer.surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.tracer_engine import *
from tracer.CoIn_rendering.rendering import *

#-------------------------Simulation parameters--------------------
# Receiver: 4 parameters cavity geometry ---------------------------------
apertureRadius = 0.34 # (m) >0
apertureDepth = 1.5 # (m) > - coneDepth
coneRadius = 0.339 # (m) >0
coneDepth = 0.01 # (m) > - apertureRadius
# Receiver: optics ------------------------------------------------
absReceiver = 0.87 # Receiver absorptivity
emsReceiver = absReceiver # Receiver emissivity
# Receiver: Temperature boundary condition -------------------------------
tempAmbiant = 20+273.15 # Ambiant temperature in K
tempReceiver = 500+273.15 # Receiver internal wall temperature in K

# Collector: paraboloidal dish geometry ----------------------------------
dishDiameter = 22. # (m) dish diameter
dishFocus = 13.1 # (m) dish focal length
# Collector: optics -----------------------------------------------
absDish = 0.06 # mirror absorptivity
sigma = 4e-3 # (rad) local x and y axis gaussian shape error distribution parameter

case = Fourparamcav(apertureRadius, apertureDepth, coneRadius, coneDepth, absReceiver, emsReceiver, dishDiameter, dishFocus, absDish, sigma)
case.VF_sim(1,1)
case_results = case.sim(Tamb=tempAmbiant, Trec=tempReceiver,CSR = 0.1, nrays=1000000 )

print  '%\n','	Spillage:', case.spill_losses,'W\n','	Blockage:', case.block_losses,'W\n','	Dish absorptivity:', case.dish_abs_losses,'W\n','	Receiver reflectivity:', case.rec_ref_losses,'W\n', 'Flux absorbed by the receiver:', case.ReceiverA,'W\n','Optical losses total:',  case.opt_losses_total,'W\n', 'thermal losses:', case.emissive_losses,'W'

#case.VF.plot()
#case_render = Renderer(case)
#case_render.show_rays()
#case.fluxmap()
#plt.show()

