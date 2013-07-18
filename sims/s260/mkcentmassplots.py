from pylab import *
import numpy

# load data
com_260 = loadtxt("s260-blobs_centerOfMass.txt")
com_261 = loadtxt("../s261/s261-blobs_centerOfMass.txt")
com_262 = loadtxt("../s262/s262-blobs_centerOfMass.txt")

# plot COM data
figure(1)
plot(com_260[:,0], com_260[:,1]-5, label='Ra = 1e4')
plot(com_261[:,0], com_261[:,1]-5, label='Ra = 1e6')
plot(com_262[:,0], com_262[:,1]-5, label='Ra = 1e8')
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
skip = 13
# compute velocity
tm_260, v_260 = calcVel(com_260[0:-1:skip,0], com_260[0:-1:skip,1])
tm_261, v_261 = calcVel(com_261[0:-1:skip,0], com_261[0:-1:skip,1])
tm_262, v_262 = calcVel(com_262[0:-1:skip,0], com_262[0:-1:skip,1])

# plot COM velocity dataZ
figure(2)
plot(tm_260, v_260, label='Ra = 1e4')
plot(tm_261, v_261, label='Ra = 1e6')
plot(tm_262, v_262, label='Ra = 1e8')
title('Center-of-Mass (radial coordinate)')
xlabel('Time [s]')
ylabel('Radial velocity')
legend(loc='upper right')
savefig('blob-vrad.png')

show()
