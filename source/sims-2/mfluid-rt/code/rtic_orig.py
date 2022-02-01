from pylab import *
import numpy

c = 3.e8
eps0 = 8.8e-12
mu0 = 4.e-7*pi
mi = 1.67e-27
q = 1.6e-19
k = 1.6e-19

n1 = 4.e31
wpi = sqrt(n1*q*q/(eps0*mi))
di = c/wpi
YL = 0.0
YU = 148*di
CELLS = 148*5
y = linspace(YL,YU,CELLS)

A = 0.538
n2 = (1-A)/(1+A)*n1
T = 1.0
beta = 0.167
b1 = sqrt(n1*k*T*2*mu0/beta)
wci = q*b1/mi
va = b1/(mu0*mi*n1)
g = 0.0115*wci*va

n = numpy.zeros([CELLS],numpy.float)
b = numpy.zeros([CELLS],numpy.float)

yc = (YU-YL)/2.0
for i in range(CELLS):
    if (y[i]<yc):
        n[i] = n1
    else:
        n[i] = n2
    if (y[i]<yc):
        b[i] = sqrt(b1**2 + 2*mu0*(n1*k*T - n[i]*k*T + mi*g*n1*(y[i]-y[0])))
    else:
        b[i] = sqrt(b1**2 + 2*mu0*(n1*k*T - n[i]*k*T + mi*g*n1*(yc-y[0]) + mi*g*n2*(y[i]-yc)))
        
figure(1), plot(y/(YU-YL),b/b[CELLS-1], label="Magnetic Field")
figure(1), plot(y/(YU-YL),n/n1, label="Density")
grid()
legend(loc="best")
show()
