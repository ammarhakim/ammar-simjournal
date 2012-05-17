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

fig = pylab.figure(2)
fig.subplots_adjust(hspace=0.4)
fig.subplots_adjust(wspace=0.4)

pylab.subplot(2,1,1)
tr_140 = pylab.loadtxt('../s140/s140-vortex-waltz_totalEnergy')
tr_141 = pylab.loadtxt('../s141/s141-vortex-waltz_totalEnergy')
tr_142 = pylab.loadtxt('../s142/s142-vortex-waltz_totalEnergy')

refTe = tr_140[0,1]

pylab.plot(tr_140[:,0], abs(tr_140[:,1]-tr_140[0,1])/tr_140[0,1], label='CFL 0.2')
pylab.plot(tr_141[:,0], abs(tr_141[:,1]-tr_141[0,1])/tr_141[0,1], label='CFL 0.1')
pylab.plot(tr_142[:,0], abs(tr_142[:,1]-tr_142[0,1])/tr_142[0,1], label='CFL 0.05')

pylab.legend(loc='upper left')
pylab.title('Energy Error History')
pylab.xlabel('Time [s]')
pylab.ylabel('Total Energy')

print "64x64", abs(tr_140[-1,1]-tr_140[0,1])
print "128x128", abs(tr_141[-1,1]-tr_141[0,1])
print "256x256", abs(tr_142[-1,1]-tr_142[0,1])

pylab.subplot(2,1,2)
tr_140 = pylab.loadtxt('../s140/s140-vortex-waltz_totalEnstrophy')
tr_141 = pylab.loadtxt('../s141/s141-vortex-waltz_totalEnstrophy')
tr_142 = pylab.loadtxt('../s142/s142-vortex-waltz_totalEnstrophy')

pylab.plot(tr_140[:,0], abs(tr_140[:,1]-tr_140[0,1])/tr_140[0,1], label='CFL 0.2')
pylab.plot(tr_141[:,0], abs(tr_141[:,1]-tr_141[0,1])/tr_141[0,1], label='CFL 0.1')
pylab.plot(tr_142[:,0], abs(tr_142[:,1]-tr_142[0,1])/tr_142[0,1], label='CFL 0.05')

print "64x64", abs(tr_140[-1,1]-tr_140[0,1])
print "128x128", abs(tr_141[-1,1]-tr_141[0,1])
print "256x256", abs(tr_142[-1,1]-tr_142[0,1])

refTe = tr_140[0,1]
pylab.legend(loc='upper left')
pylab.title('Enstrophy Error History')
pylab.xlabel('Time [s]')
pylab.ylabel('Enstrophy Error')
pylab.savefig('s140s141s142-vortex-waltz-totalEnergyEnstrophy_cmp.png')
pylab.close()

