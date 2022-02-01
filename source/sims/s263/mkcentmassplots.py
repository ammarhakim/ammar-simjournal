from pylab import *
import numpy

# load data
com_260 = loadtxt("../s260/s260-blobs_centerOfMass.txt")
com_263 = loadtxt("s263-blobs_centerOfMass.txt")

# plot COM data
figure(1)
plot(com_260[:,0], com_260[:,1]-5, label='260')
plot(com_263[:,0], com_263[:,1]-10, label='263')
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
tm_263, v_263 = calcVel(com_263[0:-1:skip,0], com_263[0:-1:skip,1])

# plot COM velocity dataZ
figure(2)
plot(tm_260, v_260, label='260')
plot(tm_263, v_263, label='263')
title('Radial velocity')
xlabel('Time [s]')
ylabel('Radial velocity')
legend(loc='upper right')
savefig('blob-vrad.png')

show()
