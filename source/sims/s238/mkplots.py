from pylab import *
import tables
import math

rc('text', usetex=True)
font = {'weight' : 'bold', 'size' : 16}
rc('font', **font)

fh = tables.openFile("s238-gemguide-5m_q_30.h5")
q = fh.root.StructGridField

mi = 1.0
me = 1.0/25.0
n0 = 1.0
B0 = 1/15.0
p0 = B0*B0*n0/12.0
pe = p0
pi = 5*p0
vte = math.sqrt(2*pe/(n0*me))
vti = math.sqrt(2*pi/(n0*mi))

nE = q[:,:,0]/me/n0
uE = q[:,:,1]/q[:,:,0]/vte
vE = q[:,:,2]/q[:,:,0]/vte
wE = q[:,:,3]/q[:,:,0]/vte

rhoI = q[:,:,5]/mi/n0
uI = q[:,:,6]/q[:,:,5]/vti
vI = q[:,:,7]/q[:,:,5]/vti
wI = q[:,:,8]/q[:,:,5]/vti

Bx = q[:,:,13]
By = q[:,:,14]

nx = ny = 768
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)
XX, YY = meshgrid(X, Y)

ylo = 4.5/25.0*nx
yup = nx-ylo+1

figure(1)
pcolormesh(X, Y[ylo:yup], transpose(nE[:,ylo:yup]))
title('Electron number density ($n/n_0$) at t$\Omega_i = 30$')
colorbar(shrink=0.675)
#streamplot(X, Y[ylo:yup], transpose(Bx[:,ylo:yup]), transpose(By[:,ylo:yup]), color='k')
axis('image')
savefig('s238-ne.png')

figure(2)
pcolormesh(X, Y[ylo:yup], transpose(vI[:,ylo:yup]))
title('Ion inflow velocity ($u_{iz}/v_{thi}$) at t$\Omega_i = 30$')
colorbar(shrink=0.675)
#streamplot(X, Y[ylo:yup], transpose(Bx[:,ylo:yup]), transpose(By[:,ylo:yup]), color='k')
axis('image')
savefig('s238-uiz.png')

figure(3)
pcolormesh(X, Y[ylo:yup], transpose(uI[:,ylo:yup]))
title('Ion outflow velocity ($u_{ix}/v_{thi}$) at t$\Omega_i = 30$')
colorbar(shrink=0.675)
#streamplot(X, Y[ylo:yup], transpose(Bx[:,ylo:yup]), transpose(By[:,ylo:yup]), color='k')
axis('image')
savefig('s238-uix.png')

figure(4)
pcolormesh(X, Y[ylo:yup], transpose(wE[:,ylo:yup]))
title('Electron out-of-plane velocity ($u_{ez}/v_{the}$) at t$\Omega_i = 30$')
colorbar(shrink=0.675)
#streamplot(X, Y[ylo:yup], transpose(Bx[:,ylo:yup]), transpose(By[:,ylo:yup]), color='k')
axis('image')
savefig('s238-uey.png')

figure(5)
pcolormesh(X, Y[ylo:yup], transpose(uE[:,ylo:yup]))
title('Electron outflow velocity ($u_{ex}/v_{the}$) at t$\Omega_i = 30$')
colorbar(shrink=0.675)
#streamplot(X, Y[ylo:yup], transpose(Bx[:,ylo:yup]), transpose(By[:,ylo:yup]), color='k')
axis('image')
savefig('s238-uex.png')

figure(6)
plot(Y, nE[nx/2,:], 'k-')
title('Electron number density ($n/n_0$) at t$\Omega_i = 30$')
axis('tight')

savefig('s238-ne-diff.png')

figure(7)
plot(Y, Bx[nx/2,:]/B0, 'k-')
plot([3,3], [-1,0.834], 'b--')
plot([3,12.5], [0.834,0.834], 'b--')
text(12.5+0.1, 0.834, '0.834', weight='bold')
title('Normalized magnetic field')
ylabel('$B(y)/B_0$')
axis('tight')
savefig('s238-bx-diff.png')

show()

