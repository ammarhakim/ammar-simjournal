import tables
import numpy
from pylab import *

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

figure(1)
plot(tm, diffByFlux)
title('Reconnected flux rate')
xlabel('T')
ylabel('d(flux)/dt')

show()
