from pylab import *

elcCharge = 1.0
elcMass = 1/25.0

def lowPass(x, dt, C):
    y = 0*x # filtered signal
    alpha = dt/(C+dt)
    y[0] = x[0]
    for i in range(1,x.shape[0]):
        y[i] = alpha*x[i] + (1-alpha)*y[i-1]
    return y

def calcDiff(T, fld):
    Td = zeros((fld.shape[0]-1,), float)
    deriv = zeros((fld.shape[0]-1,), float)
    for i in range(fld.shape[0]-1):
        deriv[i] = (fld[i+1]-fld[i])/(T[i+1]-T[i])
        Td[i] = 0.5*(T[i+1]+T[i])
    return Td, deriv

# load data
dat = loadtxt('s311-5m-gem_byFlux.txt')
T = dat[:,0]
By = dat[:,1]

dt = T[1]-T[0]
ByFilter = lowPass(By, dt, 10.0)

figure(1)
plot(T, By, '-k')

figure(2)
plot(T

show()
