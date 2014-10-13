
from tracer.surface import *
from tracer.quadric import *
from tracer.paraboloid import *
from tracer.cone import *
from tracer.cylinder import *
from tracer.flat_surface import *
from tracer.assembly import *
from tracer.optics_callables import *
from tracer.object import *
from tracer.spatial_geometry import *
from tracer.sources import *
from tracer.tracer_engine import *
from tracer.sphere_surface import *

import numpy as N
from tracer.CoIn_rendering.rendering import *

A = Assembly()
surf = AssembledObject(surfs=[Surface(SphericalRectFacet(radius = 1., lx = 0.2, ly = 0.3), LambertianReceiver(1.))], transform = None)
A.add_object(surf)
test = TracerEngine(A)

scn = Renderer(test)
scn.show_geom()
