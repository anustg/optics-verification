import numpy as N
from view_factors_3D import *
from emissive_losses import radiosity_RTVF

#example perso:
# cylinder 1m radius, 2 sections of 1m length
#areas = [N.pi, 2.*N.pi, 2.*N.pi, N.pi]
#VF = N.array([[0., 0.618, 0.210, 0.172],
#			[0.309, 0.382, 0.204, 0.105],
#			[0.105, 0.204,0.382, 0.309],
#			[0.172, 0.210, 0.618, 0.]])
#J = [459.084, 556.13, 2139, -258.8]
#cyl = Cylinder_cavity_RTVF(1.,1., [1.,1.], 5000, 0.005)
cyl = Two_N_parameters_cavity_RTVF(apertureRadius=1., frustaRadii=[1.0000001,1.], frustaDepths=[1.,1.], coneDepth=0., el_FRUs=N.array([10.,10.]), el_CON=1., num_rays=100000., precision=0.0005)

VF = cyl.VF_esperance
areas = cyl.areas
print VF
print areas
eps = N.array([1.,0.5,0.5,0.5])
Tamb = 300.
Twall = N.array([400.,500.,float('nan')])
inc_radiation = N.array([float('nan'), float('nan'),float('nan'),1000.])

AA,bb,J,E,T,q,Q = radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation)
print 'AA:', AA
print 'bb:', bb
print 'J:', J
print 'Q:', Q
print 'T:', T
print 'Q balance=', Q[0]+N.sum(Q[1:])
'''
#Holman 8th edition p470, example 8.17 cylinder
cyl = Cylinder_cavity_RTVF(0.01,0.01, [0.01, 0.01,0.01], 10000, 0.01)
#VF = N.array(([0, 0.63, 0.195, 0.075, 0.1],
#			 [0.315, 0.37, 0.2175, 0.06, 0.0375],
#			 [0.0975, 0.2175, 0.37, 0.2175, 0.0975],
#			 [0.0375, 0.06, 0.2175, 0.37, 0.315],
#			 [0.1, 0.075, 0.195, 0.63, 0]))
VF = cyl.VF
areas = cyl.areas
print VF
#areas = N.array([N.pi*(1e-2)**2,2*N.pi*1e-2*1e-2,2*N.pi*1e-2*1e-2,2*N.pi*1e-2*1e-2,N.pi*(1e-2)**2])
eps = N.array([1,0.6,0.6,0.6,0.6])
Tamb = 293
Twall = N.array([1273,1273,1273,1273])
inc_radiation = None

AA,bb,J,E,T,q,Q = radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation)

print 'Q:', Q
print 'Q balance=', Q[0]+N.sum(Q[1:])

#Suryanarayana N.V. "Engineering Heat Transfer" 1st Ed example 8.3.1 and 8.3.2
Fa = 1.-1./N.sqrt(2.)
Fb = 1.-2.*Fa
VF = N.array(([0.,Fa,Fb,Fa],[Fa,0.,Fa,Fb],[Fb,Fa,0.,Fa],[Fa,Fb,Fa,0.])) #view factors
eps = N.array([0.9, 1., 0.1, 0.8]) #emissivities
Tamb = 400.
Twall = N.array([500., 600., 0.]) # in K .If surface is known flux, set T=0
inc_radiation = N.array([N.nan, N.nan,N.nan,
5000.]) 
areas = N.array([1,1,1,1]) # in m^2. This example is per unit area
passive = [-1]

AA,bb,J,E,T,q,Q = radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation)
print 'AA', AA
print 'J', J
print 'q:', q
print 'T:', T
print N.sum(Q)

VF = N.array(([0, 0.20005, 0.20005, 0.20005, 0.20005, 0.1998],[0.20005, 0, 0.20005, 0.1998, 0.20005, 0.20005],[0.20005, 0.20005, 0, 0.20005, 0.1998, 0.20005],[0.20005, 0.1998, 0.20005, 0.20005, 0, 0.20005],[0.20005, 0.20005, 0.1998, 0.20005, 0, 0.20005],[0.1998, 0.20005, 0.20005, 0.20005, 0.20005, 0]))
Twall = N.array([20,100,100,100,100,200])
epsi=numpy.matrix('1;0.8;0.8;0.8;0.8;0.8')
Barea=0.1*numpy.matrix('1;1;1;1;1;1')
'''
