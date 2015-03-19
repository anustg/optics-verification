from math import pi, sqrt, tan
import numpy as N
import scipy.linalg as S

def radiosity_RTVF(VF, areas, eps, Tamb, Twall, inc_radiation=None):
    '''
    Solve the radiosity problem depending on the geometry and a fixed temmperature/incoming radiation boundary condition.

    Arguments:
    VF - view factor matrix of the geometry calculated
    areas - areas of the view factor elements (m2)
    eps - emissivity of the surfaces, can be an array decribing descrete emissivity values or a single number applicable to all surfaces except the aperture.
    Tamb - ambiant temperature (K)
    Twall - wall boundary condition, can be an array decribing descrete temperature values or a single number applicable to all surfaces except the aperture (K)
    inc_radiation - overrides the temperature condition if specified and specifies the ratiative power coming to the surface (W)

    Returns:
    AA - areas of the elements (m2)
    bb - incoming radiation per element (W)
    J - Radiosities (W/m2)
    E - black body emission flux (W/m2)
    T - Temepraturs of teh elements (K)
    q - Net radiative flux (W/m2)
    Q - Net radiative power (W)
    ''' 
    A = areas
    n = N.shape(VF)[0]
    sigma = 5.6677e-8 # Stefan-Boltzman constant
    if type(eps) == float:
        eps = N.hstack((1.,N.ones(n-1)*eps))
    # The radiosity problem is formulated as [AA][J]=[bb], and the following matrices are
    # defined:
    AA = N.zeros((n,n)) 
    bb = N.ones(n)

    if type(Twall)==int or type(Twall)==float:
        Twall = N.ones(n-1)*Twall

    T = N.hstack((Tamb, Twall))

    # radiosity problem, assuming specified boundary temperatures
    # The problem is formulated as [AA][J]=[bb], where [bb]=eps[i]*Eb
    for i in range(n):
        # LHS, matrix 'AA':
        for j in range(n):
            if i == j:
                AA[i,j] = 1. - VF[i,j]*(1. - eps[i])
            else:
                AA[i,j] = -VF[i,j]*(1.-eps[i])

        # RHS 'bb': for 1..n-1, setting up the temperture matrix, T[i]
        bb[i] = eps[i]*sigma*T[i]**4.
        if (inc_radiation != None):
            if ~N.isnan(inc_radiation[i]):
                bb[i] = eps[i]*inc_radiation[i]

    # Matrix inversion:
    J = N.linalg.solve(AA, bb)

    # Calculate element-wise flux density q and total flux Q.
    q = N.zeros(n)
    Q = N.zeros(n)
    E = sigma*T**4.

    for i in range(n):
        if eps[i]!=1.:
            q[i] = eps[i]/(1.-eps[i])*(E[i]-J[i]) #(W/m2)

            if (inc_radiation != None):
                if ~N.isnan(inc_radiation[i]):
                    q[i] = eps[i]/(1.-eps[i])*(inc_radiation[i]-J[i])
                    T[i] = ((A[i]*q[i]*(1.-eps[i])/(eps[i]*A[i])+J[i])/sigma)**(1./4.)

        else:
            q[i] = J[i]-N.sum(VF[i,:]*J) #(W/m2)
            if (inc_radiation != None):
                if ~N.isnan(inc_radiation[i]):
                    q[i] = (inc_radiation[i]-N.sum(VF[i,:]*J))
                    T[i] = ((A[i]*q[i]*(1.-eps[i])/(eps[i]*A[i])+J[i])/sigma)**(1./4.)  
    Q = A*q #(W)
    return AA,bb,J,E,T,q,Q
# vim: ts=4 et:
