import pylab
import tables
import math
import numpy
#pylab.rc('text', usetex=True)

dx = 1/800.
X = pylab.linspace(0+0.5*dx, 1-0.5*dx, 800)
T = pylab.linspace(0, 5e-9, 101)

for i in range(101):
    fh = tables.openFile("../s70/s70-plasmabeach_q_%d.h5" % i)
    q800 = fh.root.StructGridField
    pylab.figure(1)
    pylab.plot(X, q800[:,11], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-1.2e-11, 1.2e-11])
    pylab.title("T=%g" % T[i])
    pylab.savefig('s70-Ey_%05d.png' % i)
    pylab.close()
