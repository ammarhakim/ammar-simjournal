from pylab import *
import numpy
import tables

Lx = 100.0
Ly = 50.0
nFrame = 50
elcEnergy = numpy.zeros((nFrame+1,), numpy.float)
ionEnergy = numpy.zeros((nFrame+1,), numpy.float)
emEnergy = numpy.zeros((nFrame+1,), numpy.float)

Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s301-5m-double-periodic_q_%d.h5" % i)
    q = fh.root.StructGridField

    dx = Lx/q.shape[0]
    dy = Ly/q.shape[1]

    emEnergy[i] = 0.5*dx*dy*sum(q[:,:,13]**2 + q[:,:,14]**2 + q[:,:,15]**2)

figure(1)
plot(2*Tm, emEnergy, '-b', label='EM')
xlabel('Time')
ylabel('Magnetic energy')
savefig('s301-mag-energy.png')

show()


