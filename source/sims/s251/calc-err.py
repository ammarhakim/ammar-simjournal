from pylab import *
import tables

def exactSol(X, Y, t):
    return exp(-2*t)*sin(X)*cos(Y)

fh = tables.openFile("s251-dg-diffuse-2d_q_1.h5")
q = fh.root.StructGridField
nx, ny, nc = q.shape
dx = 2*pi/nx
Xf = linspace(0, 2*pi-dx, nx)
dy = 2*pi/ny
Yf = linspace(0, 2*pi-dy, ny)
XX, YY = meshgrid(Xf, Yf)

Xhr = linspace(0, 2*pi, 101)
Yhr = linspace(0, 2*pi, 101)
XXhr, YYhr = meshgrid(Xhr, Yhr)
fhr = exactSol(XXhr, YYhr, 1.0)

figure(1)
pcolormesh(Xhr, Yhr, fhr)
colorbar()
figure(2)
pcolormesh(Xf, Yf, q[:,:,0])
colorbar()

# compute error
fex = exactSol(XX, YY, 1.0)
error = abs(fex.transpose()-q[:,:,0]).sum()/(nx*ny);

print "%g %g" % (dx, error)

def evalSum(coeff, fields):
    res = 0.0*fields[0]
    for i in range(len(coeff)):
        res = res + coeff[i]*fields[i]
    return res

def projectOnFinerGrid_f24(Xc, Yc, q):
    dx = Xc[1]-Xc[0]
    dy = Yc[1]-Yc[0]
    nx = Xc.shape[0]
    ny = Yc.shape[0]

    # mesh coordinates
    Xn = linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 2*nx+1) # one more
    Yn = linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 2*ny+1) # one more
    XXn, YYn = meshgrid(Xn, Yn)

    # data
    qn = zeros((2*Xc.shape[0], 2*Yc.shape[0]), float)

    v1 = q[:,:,0]
    v2 = q[:,:,1]
    v3 = q[:,:,2]
    v4 = q[:,:,3]

    vList = [v1,v2,v3,v4]

    # node 1
    c1 = [0.5625,0.1875,0.0625,0.1875]
    qn[0:2*nx:2, 0:2*ny:2] = evalSum(c1, vList)

    # node 2
    c2 = [0.1875,0.5625,0.1875,0.0625]
    qn[1:2*nx:2, 0:2*ny:2] = evalSum(c2, vList)

    # node 3
    c3 = [0.1875,0.0625,0.1875,0.5625]
    qn[0:2*nx:2, 1:2*ny:2] = evalSum(c3, vList)

    # node 4
    c4 = [0.0625,0.1875,0.5625,0.1875]
    qn[1:2*nx:2, 1:2*ny:2] = evalSum(c4, vList)
   
    return XXn, YYn, qn

Xc = linspace(0.5*dx, 2*pi-0.5*dx, nx)
Yc = linspace(0.5*dy, 2*pi-0.5*dy, ny)

Xp, Yp, qp = projectOnFinerGrid_f24(Xc, Yc, q)
figure(1)
subplot(1,2,1)
pcolormesh(Xp, Yp, transpose(qp))
title('RDG t=1')
colorbar(shrink=0.5)
axis('image')
subplot(1,2,2)
pcolormesh(Xhr, Yhr, fhr)
title('Exact t=1')
colorbar(shrink=0.5)
axis('image')
savefig('s251-exact-cmp.png')

show()
