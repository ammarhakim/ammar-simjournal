import tables
import numpy
from pylab import *

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*(ionMass*n0+elcMass*n0))
ionCycl = ionCharge*B0/ionMass

def calcDeriv(T, func):
    nt = T.shape[0]-1
    tm = numpy.zeros((nt,), numpy.float)
    vx = numpy.zeros((nt,), numpy.float)

    for i in range(nt):
        tm[i] = 0.5*(T[i+1]+T[i])
        vx[i] = (func[i+1]-func[i])/(T[i+1]-T[i])

    return tm, vx

dat = loadtxt("s281-gem-tenmom_byFlux.txt")
T = dat[:,0]
byFlux = dat[:,1]
tm, diffByFlux = calcDeriv(T, byFlux)
tmNorm = tm*ionCycl

figure(1)
plot(tmNorm, diffByFlux/(LX/2.0)/B0/va)
#plot([tmNorm[0],tmNorm[-1]], [0.1,0.1], '--k')
#text(50, 0.11, '$0.1v_A$')
title('Normalized reconnection rate')
xlabel('T')
ylabel('d(flux)/dt')

show()
