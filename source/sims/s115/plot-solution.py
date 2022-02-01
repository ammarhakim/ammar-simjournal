import pylab
import tables
import math
import numpy

def projectOnFinerGrid_f(Xc, Yc, q):
    dx = Xc[1]-Xc[0]
    dy = Yc[1]-Yc[0]
    nx = Xc.shape[0]
    ny = Yc.shape[0]

    # mesh coordinates
    Xn = pylab.zeros((2*Xc.shape[0],), float)
    Xn[0:nx] = Xc-0.25*dx
    Xn[nx:] = Xc+0.25*dx
    Xn.sort()

    Yn = pylab.zeros((2*Yc.shape[0],), float)
    Yn[0:nx] = Yc-0.25*dx
    Yn[nx: ] = Yc+0.25*dx
    Yn.sort()

    qn = pylab.zeros((2*Xc.shape[0],2*Yc.shape[0]), float)

    c0 = q[:,:,0]
    c1 = q[:,:,1]
    c2 = q[:,:,2]
    c3 = q[:,:,3]
    c4 = q[:,:,4]
    c5 = q[:,:,5]
    c6 = q[:,:,6]
    c7 = q[:,:,7]

    # node 0
    qn[0:2*nx:2, 0:2*ny:2] = (9*c7)/16.0+(3*c6)/16.0+(3*c5)/16.0+(9*c4)/16.0-(3*c3)/16.0-c2/8.0-(3*c1)/16.0

    # node 1
    qn[1:2*nx:2, 0:2*ny:2] = (3*c7)/16.0+(3*c6)/16.0+(9*c5)/16.0+(9*c4)/16.0-c3/8.0-(3*c2)/16.0-(3*c0)/16.0

    # node 2
    qn[1:2*nx:2, 1:2*ny:2] = (3*c7)/16.0+(9*c6)/16.0+(9*c5)/16.0+(3*c4)/16.0-(3*c3)/16.0-(3*c1)/16.0-c0/8.0

    # node 3
    qn[0:2*nx:2, 1:2*ny:2] = (9*c7)/16.0+(9*c6)/16.0+(3*c5)/16.0+(3*c4)/16.0-(3*c2)/16.0-c1/8.0-(3*c0)/16.0

    return Xn, Yn, qn

fh = tables.openFile("s115-pb-advection-2d_chi_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# get final solution
q_1 = fh.root.StructGridField
Xn, Yn, qn_1 = projectOnFinerGrid_f(Xc, Yc, q_1)

# get intial solution
fh = tables.openFile("s115-pb-advection-2d_chi_0.h5")
q_0 = fh.root.StructGridField
Xn, Yn, qn_0 = projectOnFinerGrid_f(Xc, Yc, q_0)

nx, ny = Xn.shape[0], Yn.shape[0]

# make plot
pylab.figure(1)

pylab.subplot(1,2,1)
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1))
pylab.axis('image')

pylab.subplot(1,2,2)
pylab.plot(Xn, qn_1[:,ny/2], '-ro', Xn, qn_0[:,ny/2], '-k')

pylab.savefig('s115-projected-solution.png')

print nx, ny

# compute error
err = numpy.abs(qn_1-qn_0).sum()/(nx*ny)
#err = numpy.abs(q_1.read()-q_0.read()).sum()/(cells[0]*cells[1]*8)
print math.sqrt(dx*dy), err

pylab.show()
