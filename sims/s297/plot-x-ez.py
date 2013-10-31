from pylab import *
import numpy
import tables

Lx = 100.0
Ly = 50.0
nFrame = 30
ez = numpy.zeros((nFrame+1,), numpy.float)
jz = numpy.zeros((nFrame+1,), numpy.float)

Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s297-gem-tenmom_q_%d.h5" % i)
    q = fh.root.StructGridField
    
    nx, ny = q.shape[0], q.shape[1] 
    dx = Lx/nx
    dy = Ly/ny

    ez[i] = q[nx/2,ny/2,22]
    jz[i] = q[nx/2,ny/2,3]*25

figure(1)
plot(Tm, ez, '-b')
xlabel('Time')
ylabel('Ez')

figure(2)
plot(Tm, jz, '-b')
xlabel('Time')
ylabel('Jz')

show()


