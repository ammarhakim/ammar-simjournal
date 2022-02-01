import numpy
from pylab import *
import tables

rc('text', usetex=True)

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*elcMass*n0)
ionCycl = ionCharge*B0/ionMass

fh = tables.openFile("s282-gem-tenmom_q_30.h5")
tm = fh.root.timeData._v_attrs['vsTime']*ionCycl
q = fh.root.StructGridField
nx, ny = q.shape[0], q.shape[1]

X = linspace(-LX/2, LX/2, nx)
Y = linspace(-LY/2, LY/2, ny)
XX, YY = meshgrid(X, Y)

jz = q[:,:,3]
ux = q[:,:,1]/q[:,:,0]
ux1 = ux[:,ny/2]
minix = ux1.argmin()
maxix = ux1.argmax()
minx = X[minix]
maxx = X[maxix]
xten = maxx-minx

figure(1)
subplot(2,1,2)
pcolormesh(XX, YY, transpose(jz))
axis('tight')
plot([minx,minx], [-LY/2,LY/2], '--w')
plot([maxx,maxx], [-LY/2,LY/2], '--w')

title('Electron out-of-plane current, $t\Omega_{ci}$=%g. Extension=%g' % (tm,xten) )
#colorbar()

#figure(2)
subplot(2,1,1)
plot(X, ux1)
axis('tight')
yl = gca().get_ylim()
plot([minx,minx], yl, '--k')
plot([maxx,maxx], yl, '--k')
title('10M Electron out-flow velocity at $t\Omega_{ci}$=%g' % tm)
savefig('s282_extension.png')

figure(2)
Bfld = sqrt(q[:,:,23]**2+q[:,:,24]**2+q[:,:,25]**2)
uxNorm = ux/(Bfld/sqrt(mu0*elcMass*q[:,:,0]))
plot(X, uxNorm[:,250])
axis('tight')
title('10M Normalized electron out-flow velocity at $t\Omega_{ci}$=%g' % tm)
savefig('s282_normvx.png')

show()
