import pylab
import tables
import math
import numpy

pylab.rc('text', usetex=True)

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

fh = tables.openFile("s120-pb-advection-rb_chi_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

T = pylab.linspace(0, 4*math.pi, 41)

# for i in range(41):
#     print "Workin on %d" % i
#     fh = tables.openFile("s120-pb-advection-rb_chi_%d.h5" % i)
#     # get solution
#     q = fh.root.StructGridField
#     Xn, Yn, qn_1 = projectOnFinerGrid_f3(Xc, Yc, q)


#     pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1), vmin=0.0, vmax=0.5)
#     pylab.title('T=%f' % T[i])
#     pylab.axis('image')

#     pylab.savefig('s120-projected-solution_%05d.png' % i)
#     pylab.close()

# make subplots to put into journal

Ts = [r'0', r'$\pi/2$', r'$\pi$', r'$3\pi/2$']
count = 0
j = [1,2,3,4]
pylab.figure(1)
for i in [0,5,10,15]:
    print "Workin on %d" % i
    fh = tables.openFile("s120-pb-advection-rb_chi_%d.h5" % i)
    # get solution
    q = fh.root.StructGridField
    Xn, Yn, qn_1 = projectOnFinerGrid_f3(Xc, Yc, q)

    pylab.subplot(2,2,j[count])
    pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1), vmin=0.0, vmax=0.5)
    pylab.plot([0,1], [0.5,0.5], '-w', linewidth=0.5)
    pylab.plot([0.5,0.5], [0,1], '-w', linewidth=0.5)
    pylab.title(r'T=%s' % Ts[count])
    pylab.axis('image')
    count = count + 1

pylab.savefig('s120-snapshots.png')
pylab.close()


