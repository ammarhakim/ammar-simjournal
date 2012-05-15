import pylab
import tables
import math
import numpy

def projectOnFinerGrid(Xc, Yc, q):
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
    for i in range(nx):
        for j in range(ny):
            qn[2*i,2*j] = 1/16.0*(9*q[i,j,0]+3*q[i,j,1]+3*q[i,j,3]+q[i,j,2])

    # node 1
    for i in range(nx):
        for j in range(ny):
            qn[2*i+1,2*j] = 1/16.0*(9*q[i,j,1]+3*q[i,j,2]+3*q[i,j,0]+q[i,j,3])

    # node 2
    for i in range(nx):
        for j in range(ny):
            qn[2*i+1,2*j+1] = 1/16.0*(9*q[i,j,2]+3*q[i,j,1]+3*q[i,j,3]+q[i,j,0])

    # node 3
    for i in range(nx):
        for j in range(ny):
            qn[2*i,2*j+1] = 1/16.0*(9*q[i,j,3]+3*q[i,j,2]+3*q[i,j,0]+q[i,j,1])

    return Xn, Yn, qn

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

xlo, ylo = 0.0, 0.0
xup, yup = 2*math.pi, 2*math.pi
nx, ny = 64, 64

dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

tEnd = 8
nFrame = 11

T = pylab.linspace(0, tEnd, nFrame)

for i in range(nFrame):
    print "Working on ", i
    fh = tables.openFile("s127-double-shear_chi_%d.h5" % i)
    q = fh.root.StructGridField

    Xn, Yn, qn = projectOnFinerGrid_f3(Xc, Yc, q)

    fig = pylab.figure(1)
    pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=-5.0, vmax=5.0)
    pylab.axis('image')
    pylab.title('T=%f' % T[i])
    
    pylab.savefig('s127-double-shear_%05d.png' % i)
    pylab.close()

    fh.close()

