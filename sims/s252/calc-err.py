from pylab import *
import tables

def exactSol(X, Y, t):
    return exp(-2*t)*sin(X)*cos(Y)

fh = tables.openFile("s252-dg-diffuse-2d_q_1.h5")
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

show()
