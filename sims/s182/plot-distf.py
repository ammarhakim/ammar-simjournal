import pylab
import tables
import math

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
    Yn[0:ny] = Yc-0.25*dy
    Yn[ny: ] = Yc+0.25*dy
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

xlo, xup = -2*math.pi, 2*math.pi
ylo, yup = -6.0, 6.0
nx, ny = 32, 128

X = pylab.linspace(xlo, xup, nx+1)
Y = pylab.linspace(-6, 6, ny+1)

dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

tEnd = 20.0
nFrame = 1
T = pylab.linspace(0, tEnd, nFrame+1)
frames = [0, 1]

for i in range(nFrame+1):
    print "Working on ", i
    fh = tables.openFile("s182-landau-damping-vp_distf_%d.h5" % i)
    q = fh.root.StructGridField

    Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)

    pylab.figure(i)
    pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=0.0, vmax=0.6)
    #pylab.colorbar()
    pylab.axis('tight')
    pylab.title('T=%f' % T[i])
    fh.close()
    pylab.savefig('s182-landau-damping-vp_distf_%05d.png' % i)
    pylab.close()    

