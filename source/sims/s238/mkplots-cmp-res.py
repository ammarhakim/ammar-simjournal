from pylab import *
import tables
import math

rc('text', usetex=True)
font = {'weight' : 'bold', 'size' : 16}
rc('font', **font)

figure(1)

subplot(3,1,1)
fh = tables.openFile("../s239/s239-gemguide-5m_q_30.h5")
q = fh.root.StructGridField
nx = ny = 256
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)
XX, YY = meshgrid(X, Y)
ylo = 4.5/25.0*nx
yup = nx-ylo+1
pcolormesh(X, Y[ylo:yup], transpose(q[:,ylo:yup,3]))
axis('image')

subplot(3,1,2)
fh = tables.openFile("../s240/s240-gemguide-5m_q_30.h5")
q = fh.root.StructGridField
nx = ny = 512
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)
XX, YY = meshgrid(X, Y)
ylo = 4.5/25.0*nx
yup = nx-ylo+1
pcolormesh(X, Y[ylo:yup], transpose(q[:,ylo:yup,3]))
axis('image')

subplot(3,1,3)
fh = tables.openFile("s238-gemguide-5m_q_30.h5")
q = fh.root.StructGridField
nx = ny = 768
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)
XX, YY = meshgrid(X, Y)
ylo = 4.5/25.0*nx
yup = nx-ylo+1
pcolormesh(X, Y[ylo:yup], transpose(q[:,ylo:yup,3]))
axis('image')

savefig('s238-res-cmp-c3-30.png', dpi=300)

show()
