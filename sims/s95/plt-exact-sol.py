import pylab
import tables
import math
import numpy

sq5 = 1/math.sqrt(5.0)

def exactSol(a, b, X):
    c0 = -(1/2.0 + a/12.0 + b/30.0)
    c1 = 0.0
    return X**2/2 + a*X**4/12 + b*X**6/30 + c0*X + c1

fh = tables.openFile("s95-poisson-o4-1d_phi.h5")
q = fh.root.StructGridField
nx, nc = q.shape
Xf = pylab.linspace(0, 1, nx)
qe = q[:,0]

dx = Xf[1]-Xf[0]
Xm_1 = -sq5*0.5*dx + pylab.linspace(0.5*dx, 1-0.5*dx, nx-1)
qm_1 = q[:-1,1]

Xm_2 = sq5*0.5*dx + pylab.linspace(0.5*dx, 1-0.5*dx, nx-1)
qm_2 = q[:-1,2]

a = 2.0
b = -12.0

Xhr = pylab.linspace(0, 1, 101)
fhr = exactSol(a, b, Xhr)

# make plot comparing exact to numerical solution
pylab.plot(Xhr, fhr, '-r', Xf, qe, 'ok', Xm_1, qm_1, 'ob',
           Xm_2, qm_2, 'og')

# compute error
fex_e = exactSol(a, b, Xf)
fex_m_1 = exactSol(a, b, Xm_1)
fex_m_2 = exactSol(a, b, Xm_2)
error = (
    numpy.abs(fex_e-qe).sum()
    + numpy.abs(fex_m_1-qm_1).sum()
    + numpy.abs(fex_m_2-qm_2).sum()
    )/(nx+nx+nx-2);

print "%g %g" % (dx, error)

pylab.show()
