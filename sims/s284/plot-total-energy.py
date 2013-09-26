from pylab import *
import numpy
import tables

Lx = 2.0
Ly = 2.0
nFrame = 60
emEnergy = numpy.zeros((nFrame+1,), numpy.float)
Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s284-pulsebox-wave_q_%d.h5" % i)
    q = fh.root.StructGridField

    dx = Lx/q.shape[0]
    dy = Ly/q.shape[1]

    emEnergy[i] = 0.5*dx*dy*sum(q[:,:,0]**2 + q[:,:,1]**2 + q[:,:,2]**2 + q[:,:,3]**2 + q[:,:,4]**2 + q[:,:,5]**2)

figure(1)
plot(Tm, emEnergy, '-k', label='Total')
print "Energy lost is ", (emEnergy[-1]-emEnergy[0])/emEnergy[0]*100, " percent"

show()


