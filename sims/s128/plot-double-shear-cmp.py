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

fig = pylab.figure(1)

xlo, ylo = 0.0, 0.0
xup, yup = 2*math.pi, 2*math.pi
nx, ny = 64, 64

dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

fh = tables.openFile("../s125/s125-double-shear_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)
pylab.subplot(2, 2, 1)
pylab.title('DG2, 64x64 grid')
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=-5.0, vmax=5.0)
pylab.axis('image')

fh = tables.openFile("../s127/s127-double-shear_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f3(Xc, Yc, q)
pylab.subplot(2, 2, 3)
pylab.title('DG3, 64x64 grid')
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=-5.0, vmax=5.0)
pylab.axis('image')

nx, ny = 128, 128

dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

fh = tables.openFile("../s126/s126-double-shear_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)
pylab.subplot(2, 2, 2)
pylab.title('DG2, 128x128 grid')
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=-5.0, vmax=5.0)
pylab.axis('image')

fh = tables.openFile("s128-double-shear_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f3(Xc, Yc, q)
pylab.subplot(2, 2, 4)
pylab.title('DG3, 128x128 grid')
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=-5.0, vmax=5.0)
pylab.axis('image')

pylab.savefig('s125to128-double-shear-cmp.png')

fig = pylab.figure(2)
tr_125 = pylab.loadtxt('../s125/s125-double-shear_totalEnergy')
tr_126 = pylab.loadtxt('../s126/s126-double-shear_totalEnergy')
tr_127 = pylab.loadtxt('../s127/s127-double-shear_totalEnergy')
tr_128 = pylab.loadtxt('../s128/s128-double-shear_totalEnergy')

refTe = tr_128[0,1]
pylab.plot(tr_125[:,0], tr_125[:,1]-tr_125[0,1]+refTe, label='DG2 64x64')
pylab.plot(tr_126[:,0], tr_126[:,1]-tr_126[0,1]+refTe, label='DG2 128x128')
pylab.plot(tr_127[:,0], tr_127[:,1]-tr_127[0,1]+refTe, label='DG3 64x64')
pylab.plot(tr_128[:,0], tr_128[:,1]-tr_128[0,1]+refTe, label='DG3 128x128')
pylab.legend(loc='lower left')
pylab.title('Total Energy History')
pylab.xlabel('Time [s]')
pylab.ylabel('Total Energy')
pylab.savefig('s125to128-double-shear-totalEnergy_cmp.png')
pylab.close()

fig = pylab.figure(3)
tr_125 = pylab.loadtxt('../s125/s125-double-shear_totalEnstrophy')
tr_126 = pylab.loadtxt('../s126/s126-double-shear_totalEnstrophy')
tr_127 = pylab.loadtxt('../s127/s127-double-shear_totalEnstrophy')
tr_128 = pylab.loadtxt('../s128/s128-double-shear_totalEnstrophy')

refTe = tr_128[0,1]
pylab.plot(tr_125[:,0], tr_125[:,1]-tr_125[0,1]+refTe, label='DG2 64x64')
pylab.plot(tr_126[:,0], tr_126[:,1]-tr_126[0,1]+refTe, label='DG2 128x128')
pylab.plot(tr_127[:,0], tr_127[:,1]-tr_127[0,1]+refTe, label='DG3 64x64')
pylab.plot(tr_128[:,0], tr_128[:,1]-tr_128[0,1]+refTe, label='DG3 128x128')
pylab.legend(loc='lower left')
pylab.title('Total Enstrophy History')
pylab.xlabel('Time [s]')
pylab.ylabel('Total Enstropy')
pylab.savefig('s125to128-double-shear-totalEnstrophy_cmp.png')

fh.close()

