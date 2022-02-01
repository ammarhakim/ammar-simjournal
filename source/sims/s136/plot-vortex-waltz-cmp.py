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

tEnd = 100.0
nFrame = 11

T = pylab.linspace(0, tEnd, nFrame)

fig = pylab.figure(1)

xlo, ylo = 0.0, 0.0
xup, yup = 10, 10

nx, ny = 256, 256
dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

fh = tables.openFile("s136-vortex-waltz_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)

pylab.subplot(1,2,2)
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=0.0, vmax=1.0)
pylab.axis('image')
pylab.title('T=%f' % T[10])

nx, ny = 64, 64
dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

fh = tables.openFile("../s134/s134-vortex-waltz_chi_10.h5")
q = fh.root.StructGridField
Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)

pylab.subplot(1,2,1)
pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=0.0, vmax=1.0)
pylab.axis('image')
pylab.title('T=%f' % T[10])
    
pylab.savefig('s134s136-vortex-waltz_cmp.png')
pylab.close()

fig = pylab.figure(2)
fig.subplots_adjust(hspace=0.4)
fig.subplots_adjust(wspace=0.4)

pylab.subplot(2,1,1)
tr_134 = pylab.loadtxt('../s134/s134-vortex-waltz_totalEnergy')
tr_135 = pylab.loadtxt('../s135/s135-vortex-waltz_totalEnergy')
tr_136 = pylab.loadtxt('../s136/s136-vortex-waltz_totalEnergy')

refTe = tr_134[0,1]

pylab.plot(tr_134[:,0], abs(tr_134[:,1]-tr_134[0,1])/tr_134[0,1], label='64x64')
pylab.plot(tr_135[:,0], abs(tr_135[:,1]-tr_135[0,1])/tr_135[0,1], label='128x128')
pylab.plot(tr_136[:,0], abs(tr_136[:,1]-tr_136[0,1])/tr_136[0,1], label='256x256')
pylab.legend(loc='upper left')
pylab.title('Energy Error History')
pylab.xlabel('Time [s]')
pylab.ylabel('Total Energy')

print "64x64", abs(tr_134[-1,1]-tr_134[0,1])/tr_134[0,1]*100
print "128x128", abs(tr_135[-1,1]-tr_135[0,1])/tr_135[0,1]*100
print "256x256", abs(tr_136[-1,1]-tr_136[0,1])/tr_136[0,1]*100

pylab.subplot(2,1,2)
tr_134 = pylab.loadtxt('../s134/s134-vortex-waltz_totalEnstrophy')
tr_135 = pylab.loadtxt('../s135/s135-vortex-waltz_totalEnstrophy')
tr_136 = pylab.loadtxt('../s136/s136-vortex-waltz_totalEnstrophy')

pylab.plot(tr_134[:,0], abs(tr_134[:,1]-tr_134[0,1])/tr_134[0,1], label='64x64')
pylab.plot(tr_135[:,0], abs(tr_135[:,1]-tr_135[0,1])/tr_135[0,1], label='128x128')
pylab.plot(tr_136[:,0], abs(tr_136[:,1]-tr_136[0,1])/tr_136[0,1], label='256x256')

refTe = tr_134[0,1]
pylab.legend(loc='upper left')
pylab.title('Enstrophy Error History')
pylab.xlabel('Time [s]')
pylab.ylabel('Enstrophy Error')
pylab.savefig('s134s135s136-vortex-waltz-totalEnergyEnstrophy_cmp.png')
pylab.close()

fh.close()

