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
    Yn[0:ny] = Yc-0.25*dy
    Yn[ny: ] = Yc+0.25*dy
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

def plotLines(X, fld):
    for i in range(fld.shape[0]):
        pylab.plot([X[i], X[i+1]], [fld[i,0], fld[i,1]], '-k')

xlo, xup = -1.0, 1.0
ylo, yup = -6.0, 6.0
nx, ny = 64, 128

X = pylab.linspace(xlo, xup, nx+1)
Y = pylab.linspace(-6, 6, ny+1)

dx = (xup-xlo)/nx
dy = (yup-ylo)/ny

Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
Yc = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

tEnd = 20.0
nFrame = 100

T = pylab.linspace(0, tEnd, nFrame+1)

# for contouring
XX, VV = pylab.meshgrid(X,Y)
exv = 0.5*VV**2 + XX**2

for i in range(nFrame+1):
    print "Working on ", i
    fh = tables.openFile("s150-vlasov-fp_distf_%d.h5" % i)
    q = fh.root.StructGridField

    Xn, Yn, qn = projectOnFinerGrid_f(Xc, Yc, q)

    pylab.figure(i)
    pylab.pcolormesh(Xn, Yn, pylab.transpose(qn), vmin=0.0, vmax=0.4)
    pylab.colorbar()
    pylab.contour(XX, VV, exv, levels=[1], colors='k')
    pylab.axis('tight')
    pylab.title('T=%f' % T[i])
    pylab.savefig('s150-vlasov-fp_distf_%05d.png' % i)
    
    fh.close()

pylab.close()
