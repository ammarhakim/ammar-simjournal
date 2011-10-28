import pylab
import tables
import math
import numpy
#pylab.rc('text', usetex=True)

xlo = 0.0
xup = 0.14
nx = 200
dx = (xup-xlo)/nx
X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

driveF = 15.0e9
tEnd = 1.5e-9
nDump = 20
T = pylab.linspace(0, tEnd, nDump+1)

for i in range(nDump+1):
    fh = tables.openFile("s72-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField
    pylab.figure(1)
    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-0.1, 0.1])
    pylab.title("T=%g" % T[i])
    pylab.savefig('s70-Ey_%05d.png' % i)
    pylab.close()
