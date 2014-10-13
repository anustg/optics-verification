import numpy as N
from emissive_losses import radiosity_RTVF
'''
#Holman 8th edition p470, example 8.17 cylinder
VF = N.array(([0, 0.63, 0.195, 0.075, 0.1],
			 [0.315, 0.37, 0.2175, 0.06, 0.0375],
			 [0.0975, 0.2175, 0.37, 0.2175, 0.0975],
			 [0.0375, 0.06, 0.2175, 0.37, 0.315],
			 [0.1, 0.075, 0.195, 0.63, 0]))

print N.shape(VF)
areas = N.array([N.pi*(1e-2)**2,2*N.pi*1e-2*1e-2,2*N.pi*1e-2*1e-2,2*N.pi*1e-2*1e-2,N.pi*(1e-2)**2])
eps = N.array([1,0.6,0.6,0.6,0.6])
Tamb = 293
Twall = N.array([1273,1273,1273,1273])
inc_radiation = None
passive = None

AA,bb,J,E,T,q,Q = radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation, passive)

print 'Q:', Q
print 'Q balance=', Q[0]+N.sum(Q[1:])
'''
#Suryanarayana N.V. "Engineering Heat Transfer" 1st Ed example 8.3.1 and 8.3.2
Fa = 1.-1./N.sqrt(2.)
Fb = 1.-2.*Fa
VF = N.array(([0.,Fa,Fb,Fa],[Fa,0.,Fa,Fb],[Fb,Fa,0.,Fa],[Fa,Fb,Fa,0.])) #view factors
eps = N.array([0.9, 1., 0.1, 0.8]) #emissivities
Tamb = 400.
Twall = N.array([500., 600., 0.]) # in K .If surface is known flux, set T=0
inc_radiation = N.array([0,0,0,5000.]) 
areas = N.array([1,1,1,1]) # in m^2. This example is per unit area
passive = [-1]

AA,bb,J,E,T,q,Q = radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation, passive)
print 'AA', AA
print 'J', J
print 'q:', q
print 'T:', T
'''
VF = N.array(([0, 0.20005, 0.20005, 0.20005, 0.20005, 0.1998],[0.20005, 0, 0.20005, 0.1998, 0.20005, 0.20005],[0.20005, 0.20005, 0, 0.20005, 0.1998, 0.20005],[0.20005, 0.1998, 0.20005, 0.20005, 0, 0.20005],[0.20005, 0.20005, 0.1998, 0.20005, 0, 0.20005],[0.1998, 0.20005, 0.20005, 0.20005, 0.20005, 0]))
Twall = N.array([20,100,100,100,100,200])
epsi=numpy.matrix('1;0.8;0.8;0.8;0.8;0.8')
Barea=0.1*numpy.matrix('1;1;1;1;1;1')
'''
