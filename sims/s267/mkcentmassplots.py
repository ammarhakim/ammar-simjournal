from pylab import *
import numpy

# load data
com_266 = loadtxt("../s266/s266-blobs_centerOfMass.txt")
com_267 = loadtxt("../s267/s267-blobs_centerOfMass.txt")

# plot COM data
figure(1)
plot(com_266[:,0], com_266[:,1]-10, label='266')
plot(com_267[:,0], com_267[:,1]-10, label='267')
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
tm_266, v_266 = calcVel(com_266[0:-1:skip,0], com_266[0:-1:skip,1])
tm_267, v_267 = calcVel(com_267[0:-1:skip,0], com_267[0:-1:skip,1])

# plot COM velocity dataZ
figure(2)
plot(tm_266, v_266, label='266')
plot(tm_267, v_267, label='267')
title('Radial velocity')
xlabel('Time [s]')
ylabel('Radial velocity')
legend(loc='upper right')
savefig('blob-vrad.png')

show()
