# Generates a map showing a distribution of efficiencies from each heliostat
# within a field, in histogram form. Currently the parameters are optimised
# for the Sandia field; but could be changed easily for other fields, based
# on their dimensions. Further interpolation could be performed to generate
# a decent contour map. efficiencies.csv needs to be generated manually
# after each simulation.

import matplotlib.pyplot as plt
import numpy as N

## efficiencies.csv houses a n*3 array
lists = N.loadtxt('efficiencies.csv', delimiter=',')
lists = lists.T
x = lists[0]
y = lists[1]
weight = lists[2]

## range of the histogram, based on the number of mirrors per row
rngx = 28
rngy = 10
        
bins = [21,57]

## generate histogram
H, xbins, ybins = N.histogram2d(y, x, bins, \
    range=([-0.5,rngy], [-0.5,rngx]), weights=weight)

extent = [ybins[0], ybins[-1], xbins[-1], xbins[0]]

plt.imshow(H, extent=extent, interpolation='nearest')
plt.colorbar()
plt.clim(0.15,0.4)
plt.title('Field efficiency map')
plt.show()
