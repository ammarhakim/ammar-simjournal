from pylab import *
import numpy
import pylab
import math

def areSignDiff(a, b):
    if a>0 and b>0:
        return False
    elif a<0 and b<0:
        return False
    return True

def findextrema(fld):
    nx = fld.shape[0]
    idxLst = []
    prevSlope = fld[1]-fld[0]
    for i in range(1,nx-1):
        slope = fld[i+1]-fld[i]
        if areSignDiff(prevSlope, slope):
            idxLst.append(i)
        prevSlope = slope
    return idxLst

elcCharge = 1.0
elcMass = 1/25.0
ionCharge = 1.0
ionMass = 1.0

def lowPass(x, dt, C):
    y = 0*x # filtered signal
    alpha = dt/(C+dt)
    y[0] = x[0]
    for i in range(1,x.shape[0]):
        y[i] = alpha*x[i] + (1-alpha)*y[i-1]
    return y

def integrateE(T, u0, Ez):
    uz = 0.0*Ez
    uz[0] = u0
    for i in range(1,Ez.shape[0]-1):
        dt = T[i+1]-T[i]
        uz[i] = uz[i-1] - dt*elcCharge/elcMass*0.5*(Ez[i]+Ez[i+1])
    uz[Ez.shape[0]-1] = uz[Ez.shape[0]-2]
    return uz

def integrateI(T, u0, Ez):
    uz = 0.0*Ez
    uz[0] = u0
    for i in range(1,Ez.shape[0]-1):
        dt = T[i+1]-T[i]
        uz[i] = uz[i-1] + dt*ionCharge/ionMass*0.5*(Ez[i]+Ez[i+1])
    uz[Ez.shape[0]-1] = uz[Ez.shape[0]-2]
    return uz

def integrateUDiff(vei, T, ue0, ui0, Ez):
    vie = vei*elcMass/ionMass
    va = 0.5*(vie+vei)
    ma = (ionMass+elcMass)/(ionMass*elcMass)
    uDiff = 0.0*Ez
    uDiff[0] = ue0-ui0
    for i in range(1,Ez.shape[0]-1):
        dt = T[i+1]-T[i]
        Eh = 0.5*(Ez[i+1]+Ez[i])
        uDiff[i] = 1/(1+dt*va)*(uDiff[i-1]*(1-dt*va) - dt*elcCharge*Eh*ma)

    uDiff[Ez.shape[0]-1] = uDiff[Ez.shape[0]-2]
    return uDiff

def diffU(T, u, mass, charge):
    eDiff = zeros((T.shape[0]-1,), float)
    TDiff = 0.0*eDiff
    for i in range(0,T.shape[0]-1):
        eDiff[i] = mass/charge*(u[i+1]-u[i])/(T[i+1]-T[i])
        TDiff[i] = 0.5*(T[i+1]+T[i])
    return TDiff, eDiff

# load data
dat = loadtxt('s320-gem-tenmom_xpointEz.txt')
T = dat[:,0]
Ez = dat[:,1]

dat = loadtxt('s320-gem-tenmom_xpointUze.txt')
zMomE = dat[:,1]

dat = loadtxt('s320-gem-tenmom_xpointNe.txt')
rhoE = dat[:,1]
# assume quasi-neutrality
rhoI = rhoE/elcMass*ionMass

uze = zMomE/rhoE

dt = T[1]-T[0]
EzFilter = lowPass(Ez, dt, 4.0)

figure(6)
# compute an estimate of numerical diffusion
TDiff, EzDiff = diffU(T, uze, elcMass, -elcCharge)
plot(TDiff[:-1000], EzDiff[:-1000], 'k-', label='Inferred')
plot(T[:-1000], Ez[:-1000], 'r-', label='Measured')
legend(loc='best')
title('EzDiff and Ez')

show()
