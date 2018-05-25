from math import pi, sqrt, tan
import numpy as N
import scipy.linalg as S

def radiosity_RTVF(VF, areas, eps, T= None, inc_radiation=None, q_net=None):
    '''
    Solve the radiosity problem depending on the geometry and a specified boundary condition.

    Arguments:
    VF - view factor matrix of the geometry calculated
    areas - areas of the view factor elements (m2)
    eps - emissivity of the surfaces, can be an array decribing descrete emissivity values or a single number applicable to all surfaces except the aperture.
    Tamb - ambiant temperature (K)
    Twall - wall boundary condition, can be an array decribing descrete temperature values or a single number applicable to all surfaces except the aperture (K)
    inc_radiation - specifies the ratiative power coming to the surface (W/m2) and is used if no temperature is given to determien the thermal emissions.
    q_net - specifies the net radiative heat enforced in the surface. q_net is subtracted from inc_radiation if the latter is declared.

    Returns:
    J - Radiosities (W/m2)
    E - black body emission flux (W/m2)
    T - Temperatures of the elements (K)
    q - Net radiative flux (W/m2) fraction in the surface after absorption and emission.
    Q - Net radiative power (W)
    ''' 
    A = areas
    n = N.shape(VF)[0]
    sigma = 5.6677e-8 # Stefan-Boltzman constant

    if len(eps) != len(areas):
        raise AttributeError
    if (T ==None) and (inc_radiation == None):
        raise AttributeError
    # The radiosity problem is formulated as [AA][J]=[bb], and the following matrices are
    # defined:
    AA = N.zeros((n,n)) 
    bb = N.zeros(n)

    # radiosity problem, assuming specified boundary temperatures
    # The problem is formulated as [AA][J]=[bb], where [bb]=eps[i]*Eb
    for i in range(n):

        if (inc_radiation != None):
        # LHS, matrix 'AA':
            if ~N.isnan(inc_radiation[i]):      
                for j in range(n):
                    if i == j:
                        AA[i,j] = 1. - VF[i,j]#*(1. - eps[i])
                    else:
                        AA[i,j] = - VF[i,j]#*(1.-eps[i])
                bb[i] = inc_radiation[i]
                if q_net != None:
                    if ~N.isnan(q_net[i]):
                        bb[i] = bb[i]-q_net[i]

        if (T != None):
            if ~N.isnan(T[i]):
            # LHS, matrix 'AA':
                for j in range(n):
                    if i == j:
                        AA[i,j] = 1. - VF[i,j]*(1. - eps[i])
                    else:
                        AA[i,j] = - VF[i,j]*(1.-eps[i])
                # RHS 'bb': for 1..n-1,
                bb[i] = eps[i]*sigma*T[i]**4.

    if (N.isnan(bb).any()):
        raise AttributeError('Wrong left hand side')

    # Matrix inversion:
    J = N.linalg.solve(AA, bb)

    # Calculate element-wise flux density q and total flux Q.
    q = N.zeros(n)
    Q = N.zeros(n)
    E = N.zeros(n)

    for i in range(n):
        if T != None:
            if ~N.isnan(T[i]):
                E[i] = sigma*T[i]**4.
                if eps[i] !=1.:
                    q[i] = eps[i]/(1.-eps[i])*(E[i]-J[i]) #(W/m2)
                else:
                    q[i] = eps[i]*(-E[i]+N.sum(VF[:,i]*J)) # (W/m2)
            
        if (inc_radiation != None):
            if ~N.isnan(inc_radiation[i]):
                T[i] = ((q[i]*(1.-eps[i])/eps[i]+J[i])/sigma)**0.25
                q[i] = -J[i]+eps[i]*N.sum(VF[:,i]*J)+inc_radiation[i] #(W/m2)

    E = sigma*T**4.
    Q = A*q #(W)

    return AA,bb,J,E,T,q,Q
# vim: ts=4 et:
