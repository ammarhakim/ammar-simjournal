from pylab import *
import numpy

skip = 1
# load data
comx = loadtxt("../s264/s264-blobs_centerOfMass.txt")
comy = loadtxt("../s264/s264-blobs_centerOfMass_yloc.txt")

def calcVel(T, xloc):
    nt = T.shape[0]-1
    tm = numpy.zeros((nt,), numpy.float)
    vx = numpy.zeros((nt,), numpy.float)

    for i in range(nt):
        tm[i] = 0.5*(T[i+1]+T[i])
        vx[i] = (xloc[i+1]-xloc[i])/(T[i+1]-T[i])

    return tm, vx

tmx, ux = calcVel(comx[:,0], comx[:,1])
tmy, uy = calcVel(comy[:,0], comy[:,1])

def integrateForward(x0, tm, vel):
    xnew = numpy.zeros((tm.shape[0]+1,), numpy.float)
    xnew[0] = x0
    for i in range(1,vel.shape[0]):
        xnew[i] = xnew[i-1]+vel[i-1]*(tm[i]-tm[i-1])

    return xnew

# compute total velocity at each point
vel = sqrt(ux**2+uy**2)
xnew = integrateForward(0.0, tmx, vel)

figure(1)
plot(comx[:,0], comx[:,1]-10, 'k-')
plot(comx[:,0], xnew, 'r-')

figure(2)
plot(tmx, ux, '-k')
tmm, uxx = calcVel(comx[0:-2:13,0], xnew[0:-2:13])
plot(tmm, uxx, '-r')

show()
