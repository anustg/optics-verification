import os
import numpy as N
from run import *

'''
# A1
for i in xrange(2):
    for j in xrange(3):
        case='A1.%s.%s'%(i+1,j+1)
        rt=RayTracing(num_rays=10000, rendering=False, case=case)
        rt.trace()	

# A2.1 and A2.2
for i in xrange(2):
    case='A2.%s'%(i+1)
    rt=RayTracing(num_rays=10000, rendering=False, case=case)
    rt.trace()	



# A2.3
for i in xrange(3):
    case='A2.3.%s'%(i+1)
    rt=RayTracing(num_rays=10000, rendering=False, case=case)
    rt.trace()	



# A3.1 and A3.2
for i in xrange(2):
    case='A3.%s'%(i+1)
    rt=RayTracing(num_rays=10000, rendering=False, case=case)
    rt.trace()	


# Round B

for h in xrange(2):
    for i in xrange(2):
        for j in xrange(4):
            case='B%s.%s.%s'%(h+1,i+1,j+1)
            rt=RayTracing(num_rays=10000, rendering=False, case=case)
            rt.trace()	

'''
for i in xrange(2):
    for j in xrange(1,2):
        case='C%s.%s'%(i+1,j+1)
        rt=RayTracing(num_rays=10000, rendering=False, case=case)
        rt.trace()	



