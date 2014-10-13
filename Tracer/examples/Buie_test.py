from tracer.sources import *
import matplotlib.pyplot as plt
import numpy as N

sourceCenter = N.array([[0,0,2.]]).T # Source center position
sourceDirection = N.array([0,0,-1.]) # Source normal direction
sourceRadius = 0.6 # m, Source radius
sourceAngle = 4.65e-3 # radians, sun rays angular range
G = 1000. # W/m2, solar flux

for i in N.arange(0.1,0.8,0.1):
	buie = buie_sunshape(100000, sourceCenter, sourceDirection, sourceRadius, i, G)
plt.show()
