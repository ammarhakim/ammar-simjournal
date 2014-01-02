from pylab import *
import math

# normalized beta
beta = 10.0
kpar = 0.5

vA = 1/math.sqrt(beta)
dat = loadtxt('b-10-kperp-scan')
kperp = dat[:,0]
Tp = dat[:,1]
wp = 0.5*(2*pi/Tp)/(kpar*vA)
plot(kperp, wp, 'ro')

kperpL = linspace(0.01, 1.0, 100)
wpNormEx = 1/sqrt(1+kperpL**2*vA**2)
plot(kperpL, wpNormEx, '-b')

show()
