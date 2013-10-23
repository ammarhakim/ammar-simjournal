from pylab import *
import numpy
import tables

Lx = 8*pi
Ly = 4*pi
nFrame = 50
ez = numpy.zeros((nFrame+1,), numpy.float)
jz = numpy.zeros((nFrame+1,), numpy.float)

Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s311-5m-gem_q_%d.h5" % i)
    q = fh.root.StructGridField
    
    nx, ny = q.shape[0], q.shape[1] 
    dx = Lx/nx
    dy = Ly/ny

    ez[i] = q[nx/2,ny/2,12]
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


