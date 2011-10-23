import pylab
import tables
import math
import numpy
pylab.rc('text', usetex=True)

NT = 100
Ey_tx = numpy.zeros( (NT+1, 400), numpy.float )
for i in range(NT+1):
    fh = tables.openFile("s69-plasmabeach_q_%d.h5" % i)
    q = fh.root.StructGridField
    Ey = q[:,11]
    Ey_tx[i,:] = Ey

dx = 1/400.0
X = pylab.linspace(0.5*dx, 1-0.5*dx, 400)
T = pylab.linspace(0, 5e-9, NT+1)
TT, XX = pylab.meshgrid(T, X)

# compute cutoff location
dx100 = 1/100.
deltaT = dx100/2.99792458e8
driveOmega = 3.14159265358979323846264338328/10/deltaT
xcutoff = 1-math.pow(driveOmega/25*deltaT, 1/5.0)

pylab.pcolormesh(TT, XX, Ey_tx.transpose())
pylab.plot([0, 5e-9], [xcutoff, xcutoff], 'k--', linewidth=2)
pylab.xlabel('Time [s]')
pylab.ylabel('X [m]')
pylab.title(r'$E_y(t,x)$')
pylab.savefig('s69-plasmabeach_Ey.png')
pylab.show()
