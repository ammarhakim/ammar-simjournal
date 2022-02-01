from pylab import *
import numpy

# load data
com_264 = loadtxt("../s264/s264-blobs_centerOfMass.txt")
arcOn_tm = loadtxt('craig-time.txt')
arcOn_xl = loadtxt('craig-com.txt')
arcOn_vl = loadtxt('craig-com-vel.txt')

# plot COM data
figure(1)
plot(com_264[:,0], com_264[:,1]-10, label='Gkeyll')
plot(arcOn_tm, arcOn_xl, label='arcOn')

title('Center-of-Mass (radial coordinate)')
xlabel('Time [s]')
ylabel('X location')
legend(loc='upper left')
savefig('arcOn-gkeyll-blob-xloc.png')

def calcVel(T, xloc):
    nt = T.shape[0]-1
    tm = numpy.zeros((nt,), numpy.float)
    vx = numpy.zeros((nt,), numpy.float)

    for i in range(nt):
        tm[i] = 0.5*(T[i+1]+T[i])
        vx[i] = (xloc[i+1]-xloc[i])/(T[i+1]-T[i])

    return tm, vx

# this skip is needed to smooth out some of the spiky data
skip = 1
# compute velocity
tm_264, v_264 = calcVel(com_264[0:-1:skip,0], com_264[0:-1:skip,1])

# plot COM velocity dataZ
figure(2)
plot(tm_264, v_264, label='Gkeyll')
plot(arcOn_tm[:-1], arcOn_vl, label='arcOn')
title('Radial velocity')
xlabel('Time [s]')
ylabel('Radial velocity')
legend(loc='upper right')
savefig('arcOn-gkeyll-blob-vrad.png')

show()
