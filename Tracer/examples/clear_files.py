# Resets the data files for x, y coordinates of hits at the receiver, as well
# as the energies of each hit, and cumulative hits on heliostats. Currently
# optimised for the Sandia field with 218 heliostats.

import numpy as N

N.save('x_list.npy', [])
N.save('y_list.npy', [])
N.save('energy_list.npy', [])
N.save('hits_list.npy', [0.]*218)
N.save('total_list.npy', [0.]*218)
