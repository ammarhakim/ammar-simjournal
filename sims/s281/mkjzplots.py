import numpy
from pylab import *
import tables

rc('text', usetex=True)

def findXloc(X, fx, val):
    cross = []
    for i in range(X.shape[0]-1):
        if val>fx[i] and val<fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))
        elif val<fx[i] and val>fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))            
    return cross

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

fh = tables.openFile("s281-gem-tenmom_q_30.h5")
tm = fh.root.timeData._v_attrs['vsTime']*ionCycl
q = fh.root.StructGridField
nx, ny = q.shape[0], q.shape[1]

X = linspace(-LX/2, LX/2, nx)
Y = linspace(-LY/2, LY/2, ny)
XX, YY = meshgrid(X, Y)

jz = q[:,:,3]
ux = q[:,:,1]/q[:,:,0]
ux1 = ux[:,ny/2]
uy = q[:,:,2]/q[:,:,0]
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
savefig('s281_extension.png')

figure(2)
Bfld = sqrt(q[:,:,23]**2+q[:,:,24]**2+q[:,:,25]**2)
uxNorm = ux/(Bfld/sqrt(mu0*elcMass*q[:,:,0]))
plot(X, uxNorm[:,ny/2])
axis('tight')
title('10M Normalized electron out-flow velocity at $t\Omega_{ci}$=%g' % tm)
savefig('s281_normvx.png')

figure(3)
Bfld = sqrt(q[:,:,23]**2+q[:,:,24]**2+q[:,:,25]**2)
uxNorm = ux/(Bfld/sqrt(mu0*elcMass*q[:,:,0]))
pcolormesh(XX, YY, transpose(uxNorm))
colorbar()
axis('image')
title('10M Normalized electron out-flow velocity at $t\Omega_{ci}$=%g' % tm)
savefig('s281_normvx_pcolor.png')

figure(4)
plot(Y, uy[nx/2,:])
jzHalfMax = jz[jz[:,ny/2].argmax(), ny/2]/2.0
cross = findXloc(Y, jz[nx/2,:], jzHalfMax)
plot([cross[0],cross[0]], gca().get_ylim(), '--k')
plot([cross[1],cross[1]], gca().get_ylim(), '--k')
gca().set_xlim(-1.0, 1.0)
#axis('tight')
title('10M in-flow velocity at $t\Omega_{ci}$=%g' % tm)
xlabel('Y (width %g)' % (cross[1]-cross[0]))
ylabel('Inflow velocity')
savefig('s281_vy.png')

figure(5)
uyNorm = uy/(Bfld/sqrt(mu0*elcMass*q[:,:,0]))
pcolormesh(XX, YY, transpose(uy))
colorbar()
axis('image')
title('10M in-flow velocity at $t\Omega_{ci}$=%g' % tm)
xlabel('Y')
ylabel('Inflow velocity')
savefig('s281_vy_pcolor.png')

show()
