from pylab import *
import numpy

Lx = 100.0
Ly = 50.0
NX = 200
NY = 100

B0 = 0.1
me = 1.0
mi = me*25.0
qe = -1.0
qi = 1.0
dlambda = 1.0
n0 = 1.0
ninf = 0.2*n0
psi0 = B0

dx = Lx/NX
dy = Ly/NY
X = linspace(0.5*dx, Lx-0.5*dx, NX)
Y = linspace(0.5*dy, Ly-0.5*dy, NY)
XX, YY = meshgrid(X, Y)

Bx = numpy.zeros((NX, NY), numpy.float)
n = numpy.zeros((NX, NY), numpy.float)
dBx1 = numpy.zeros((NX, NY), numpy.float)
dBy1 = numpy.zeros((NX, NY), numpy.float)
dBx2 = numpy.zeros((NX, NY), numpy.float)
dBy2 = numpy.zeros((NX, NY), numpy.float)

for i in range(NX):
    for j in range(NY):
        Bx[i,j] = B0*(-1+tanh((Y[j]-Ly/4)/dlambda)-tanh((Y[j]-3*Ly/4)/dlambda))
        n[i,j] = n0/cosh((Y[j]-Ly/4)/dlambda)**2+n0/cosh((Y[j]-3*Ly/4)/dlambda)**2+ninf

        dBx1[i,j] = -psi0*(pi/Ly)*cos(2*pi*(X[i]-Lx/4)/Lx)*sin(pi*(Y[j]-Ly/4)/Ly)
        dBy1[i,j] = psi0*(2*pi/Lx)*sin(2*pi*(X[i]-Lx/4)/Lx)*cos(pi*(Y[j]-Ly/4)/Ly)

        dBx2[i,j] = -psi0*(pi/Ly)*cos(2*pi*(X[i]+Lx/4)/Lx)*sin(pi*(Y[j]+Ly/4)/Ly)
        dBy2[i,j] = psi0*(2*pi/Lx)*sin(2*pi*(X[i]+Lx/4)/Lx)*cos(pi*(Y[j]+Ly/4)/Ly)

figure(1)
pcolormesh(XX, YY, transpose(Bx))
title('Bx(x,y)')
colorbar()

figure(2)
pcolormesh(XX, YY, transpose(n))
title('n(x,y)')
colorbar()

figure(3)
plot(Y, Bx[NX/2,:], 'r-')
xlabel('Y')
ylabel('Bx')
title('Bx(y)')

figure(4)
plot(Y, n[NX/2,:], 'r-')
xlabel('Y')
ylabel('n')
title('n(y)')

figure(7)
Bxt = Bx+dBx1+dBx2
Byt = dBy1+dBy2
Btot = sqrt(Bxt**2+Byt**2)
#contour(XX, YY, transpose(Btot))
streamplot(X, Y, transpose(Bxt), transpose(Byt), density=2)

show()

