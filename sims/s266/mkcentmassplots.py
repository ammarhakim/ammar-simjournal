from pylab import *
import numpy

# load data
com_264 = loadtxt("../s264/s264-blobs_centerOfMass.txt")
com_266 = loadtxt("../s266/s266-blobs_centerOfMass.txt")

# plot COM data
figure(1)
plot(com_264[:,0], com_264[:,1]-10, label='264')
plot(com_266[:,0], com_266[:,1]-10, label='266')
title('Center-of-Mass (radial coordinate)')
xlabel('Time [s]')
ylabel('X location')
legend(loc='upper left')
savefig('blob-xloc.png')

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
tm_266, v_266 = calcVel(com_266[0:-1:skip,0], com_266[0:-1:skip,1])

# plot COM velocity dataZ
figure(2)
plot(tm_264, v_264, label='264')
plot(tm_266, v_266, label='266')
title('Radial velocity')
xlabel('Time [s]')
ylabel('Radial velocity')
legend(loc='upper right')
savefig('blob-vrad.png')

show()
