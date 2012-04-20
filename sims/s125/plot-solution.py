import pylab
import tables
import math
import numpy

#pylab.rc('text', usetex=True)
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

def projectOnFinerGrid_f3(Xc, Yc, q):
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

fh = tables.openFile("s125-pb-advection-sf_chi_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# for comparison
# get intial solution
fh = tables.openFile("s125-pb-advection-sf_chi_0.h5")
q = fh.root.StructGridField
Xn, Yn, qn_0 = projectOnFinerGrid_f3(Xc, Yc, q)

# get final solution
fh = tables.openFile("s125-pb-advection-sf_chi_1.h5")
q = fh.root.StructGridField
Xn, Yn, qn_1 = projectOnFinerGrid_f3(Xc, Yc, q)

pylab.figure(0)
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1))
pylab.colorbar()
pylab.axis('image')
pylab.title('T=3.0')
pylab.savefig('s125-solution.png')

nx, ny = Xn.shape[0], Yn.shape[0]

# make plot
pylab.figure(1)

pylab.plot(Xn, qn_1[:,ny/2], '-ro', Xn, qn_0[:,ny/2], '-k')

pylab.savefig('s125-projected-solution.png')

# compute error
err = numpy.abs(qn_1-qn_0).sum()/(nx*ny)
print math.sqrt(dx*dy), err

# plot enstrophy
pylab.figure(2)
ensData = pylab.loadtxt("totalEnstrophy")
pylab.plot(ensData[:,0], ensData[:,1], 'r-')
pylab.xlabel('Time [s]')
pylab.ylabel('Total Enstrophy')
pylab.savefig('s125-total-enstrophy.png')

pylab.close()


