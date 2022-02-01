from pylab import *

fEnergy = loadtxt("m2-two-stream_fieldEnergy.dat")
pEnergy = loadtxt("m2-two-stream_totalPtclEnergy.dat")

T = fEnergy[:,0]
tEnergy = fEnergy[:,1] + pEnergy[:,1]
plot(T, tEnergy, '-r')
tMax = tEnergy.max()
gca().set_ylim([0.99*tMax, 1.01*tMax])
error = 100*(tEnergy[-1]-tEnergy[0])/tEnergy[0]
title('Error %g %s' % (error, "%"))
show()
