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

# load data
dat = loadtxt('s311-5m-gem_xpointEz.txt')
T = dat[:,0]
Ez = dat[:,1]

dat = loadtxt('s311-5m-gem_xpointUze.txt')
zMomE = dat[:,1]

dat = loadtxt('s311-5m-gem_xpointNe.txt')
rhoE = dat[:,1]

uz = zMomE/rhoE

dt = T[1]-T[0]
EzFilter = lowPass(Ez, dt, 3.0)

figure(1)
plot(T, Ez, '-k')
plot(T, EzFilter, 'r-', linewidth=2)
title('Ez at X-point')
xlabel('Time')
ylabel('Electric Field')
savefig('s311-ez-filtered.png')

uzInt = integrateE(T, uz[0], Ez)
figure(2)
plot(T, uzInt, '-k', label='Integrated')
plot(T, uz, '-r', label='Measured')
legend(loc='best')
title('Out of plane electron velocity at X-point')
xlabel('Time')
ylabel('Uze')
savefig('s311-uz-measured-vs-integrated.png')

show()
