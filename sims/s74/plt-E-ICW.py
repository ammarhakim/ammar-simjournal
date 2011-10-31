import pylab
import tables
import math
import numpy
#pylab.rc('text', usetex=True)

xlo = 2.83
xup = 3.17
nx = 400
dx = (xup-xlo)/nx
X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

for i in range(101):
    fh = tables.openFile("s74-icw_EM_%d.h5" % i)
    q = fh.root.StructGridField

    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.5)

    pylab.subplot(3, 1, 1)
    pylab.plot(X, q[:,0], 'k-')
    pylab.ylabel('Ex')
    currAx = pylab.gca()
    pylab.axis('tight')
    currAx.set_ylim([-0.8, 0.8])

    pylab.subplot(3, 1, 2)
    pylab.plot(X, q[:,1], 'k-')
    pylab.ylabel('Ey')
    currAx = pylab.gca()
    pylab.axis('tight')
    currAx.set_ylim([-0.8, 0.8])

    pylab.subplot(3, 1, 3)
    pylab.plot(X, q[:,2], 'k-')
    pylab.ylabel('Ez')
    currAx = pylab.gca()
    pylab.axis('tight')
    currAx.set_ylim([-0.8, 0.8])
    
    pylab.savefig('s74-ICW-E_%05d.png' % i)
    pylab.close()

