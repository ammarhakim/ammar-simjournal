import pylab
import tables
import math
import numpy
pylab.rc('text', usetex=True)

xlo = 0.0
xup = 0.14
nx = 200
dx = (xup-xlo)/nx
X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

driveF = 15.0e9
tEnd = 1.5e-9
nDump = 20
T = pylab.linspace(0, tEnd, nDump+1)

fig = pylab.figure(1)
fig.subplots_adjust(hspace=0.5)

for i in [1, 2, 3, 4]:
    fh = tables.openFile("s72-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField
    pylab.subplot(4, 1, i)
    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-0.05, 0.05])
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_x$')

pylab.savefig('s72-Ex.png')

fig = pylab.figure(2)
fig.subplots_adjust(hspace=0.5)

for i in [1, 2, 3, 4]:
    fh = tables.openFile("s72-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField
    pylab.subplot(4, 1, i)
    pylab.plot(X, q[:,11], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-1.0, 1.0])
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_y$')    

pylab.savefig('s72-Ey.png')
pylab.show()
pylab.close()
