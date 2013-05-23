from pylab import *
import tables
import math

rc('text', usetex=True)
font = {'weight' : 'bold', 'size' : 16}
rc('font', **font)

nx = ny = 768
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)
XX, YY = meshgrid(X, Y)

ylo = 4.5/25.0*nx
yup = nx-ylo+1

fr = [30,45]
tm = [30,45]

count = 0
for f in fr:
    fh = tables.openFile("s238-gemguide-5m_q_%d.h5" % f)
    elcCurr = fh.root.StructGridField[:,:,3]

    pcolormesh(X, Y[ylo:yup], transpose(elcCurr[:,ylo:yup]))
    colorbar(shrink=0.675)
    title('Electron out-of-plane current at t$\Omega_i = %g$' % tm[count])    
    axis('image')
    savefig('s238-elcCurr_%d.png' % f)
    count = count + 1
    close()


fh = tables.openFile("s238-gemguide-5m_q_45.h5")
elcCurr = fh.root.StructGridField[:,:,3]
xlo = nx/2+7.0/25.0*nx
xup = nx

ylo = ny/2-1.0/25.0*ny
yup = ny/2+1.0/25.0*ny
pcolormesh(X[xlo:xup], Y[ylo:yup], transpose(elcCurr[xlo:xup,ylo:yup]))
title('Zoom of electron out-of-plane current at t$\Omega_i = 45$')
axis('image')
savefig('s238-elcCurr_zoom_45.png')
show()
