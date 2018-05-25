# Generates a post-processed flux map from raw data over n iterations by
# dividing the energy at the receiver in each bin by the number of iterations.

import matplotlib.pyplot as plt
import numpy as N

iters = input('Number of iterations ')

rngx = 5.5
rngy = 5.5
        
bins = 50

H, xbins, ybins = N.histogram2d(N.load('x_list.npy'), N.load('y_list.npy'), bins, \
    range=([-rngx,rngx], [-rngy,rngy]), weights=N.load('energy_list.npy')/iters)

print('H', H, 'length', len(H))
print('xbins', xbins, 'length', len(xbins))
print('ybins', ybins, 'length', len(ybins))
print('weights', N.load('energy_list.npy'), 'length', len(N.load('energy_list.npy')))
extent = [ybins[0], ybins[-1], xbins[-1], xbins[0]]
plt.imshow(H, extent=extent, interpolation='nearest')
plt.colorbar()
plt.title('Flux Map Title')
plt.show()
