from pylab import *
import numpy
import tables

Lx = 8*pi
Ly = 4*pi
nFrame = 40
elcEnergy = numpy.zeros((nFrame+1,), numpy.float)
ionEnergy = numpy.zeros((nFrame+1,), numpy.float)
emEnergy = numpy.zeros((nFrame+1,), numpy.float)

Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s304-5m-gem_q_%d.h5" % i)
    q = fh.root.StructGridField
    
    nx, ny = q.shape[0], q.shape[1] 
    dx = Lx/nx
    dy = Ly/ny

    emEnergy[i] = q[nx/2,ny/2,12]

figure(1)
plot(Tm*10, emEnergy, '-b', label='EM')
xlabel('Time')
ylabel('Ez')
savefig('s304-ez.png')

show()


