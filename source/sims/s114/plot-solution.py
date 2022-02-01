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

    # node 0
    qn[0:2*nx:2, 0:2*ny:2] = 1/16.0*(9*q[:,:,0]+3*q[:,:,1]+3*q[:,:,3]+q[:,:,2])

    # node 1
    qn[1:2*nx:2, 0:2*ny:2] = 1/16.0*(9*q[:,:,1]+3*q[:,:,2]+3*q[:,:,0]+q[:,:,3])

    # node 2
    qn[1:2*nx:2, 1:2*ny:2] = 1/16.0*(9*q[:,:,2]+3*q[:,:,1]+3*q[:,:,3]+q[:,:,0])

    # node 3
    qn[0:2*nx:2, 1:2*ny:2] = 1/16.0*(9*q[:,:,3]+3*q[:,:,2]+3*q[:,:,0]+q[:,:,1])

    return Xn, Yn, qn

fh = tables.openFile("s114-pb-advection-2d_chi_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# get final solution
q = fh.root.StructGridField
Xn, Yn, qn_1 = projectOnFinerGrid_f(Xc, Yc, q)

# get intial solution
fh = tables.openFile("s114-pb-advection-2d_chi_0.h5")
q = fh.root.StructGridField
Xn, Yn, qn_0 = projectOnFinerGrid_f(Xc, Yc, q)

nx, ny = Xn.shape[0], Yn.shape[0]

# make plot
pylab.figure(1)

pylab.subplot(1,2,1)
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1))
pylab.axis('image')

pylab.subplot(1,2,2)
pylab.plot(Xn, qn_1[:,ny/2], '-ro', Xn, qn_0[:,ny/2], '-k')

pylab.savefig('s114-projected-solution.png')

# compute error
err = numpy.abs(qn_1-qn_0).sum()/(nx*ny)
print math.sqrt(dx*dy), err

pylab.show()
