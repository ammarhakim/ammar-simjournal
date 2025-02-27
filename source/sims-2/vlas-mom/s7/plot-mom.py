import gkedgbasis
import gkedata
import gkedginterpdat
from pylab import *

import pylab
import tables
import math
import numpy
import pylab
import numpy
from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 18
rcParams['xtick.major.size']           = 8 # default is 4
rcParams['xtick.major.width']          = 3 # default is 0.5
rcParams['ytick.major.size']           = 8 # default is 4
rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['xtick.minor.size']           = 4 # default is 4
rcParams['xtick.minor.width']          = 0.5 # default is 0.5
rcParams['ytick.minor.size']           = 4 # default is 4
rcParams['ytick.minor.width']          = 0.5 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
#rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

class Counter:
    def __init__(self):
        self.count = 0

    def bump(self):
        self.count = self.count+1
        return self.count
        
cnt = Counter()

def getXv(Xc, Vc):
    dx = (Xc[0,-1]-Xc[0,0])/(Xc.shape[1]-1)
    dv = (Vc[-1,0]-Vc[0,0])/(Vc.shape[0]-1)

    X1 = linspace(Xc[0,0]+0.5*dx, Xc[0,-1]-0.5*dx, Xc.shape[1]-1)
    V1 = linspace(Vc[0,0]+0.5*dv, Vc[-1,0]-0.5*dv, Vc.shape[0]-1)

    return X1, V1

# density
d = gkedata.GkeData("s7-2x3v-mom_numDensity.h5")
dg1Num = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Vc, num = dg1Num.project(0)
Xv, Vv = getXv(Xc, Vc)

Xhr = linspace(Xv[0], Xv[-1], 200) # for plotting
n = sin(2*pi*Xhr)
ux = 0.1*cos(2*pi*Xhr)
uy = 0.2*sin(2*pi*Xhr)
uz = 0.1*cos(2*pi*Xhr)
Txx = 0.75 + 0.25*cos(2*pi*Xhr)
Tyy = 0.75 + 0.25*sin(2*pi*Xhr)
Tzz = 0.75 + 0.1*sin(2*pi*Xhr)
Txy = 0.5 + 0.1*sin(2*pi*Xhr)
Txz = 0.25 + 0.1*sin(2*pi*Xhr)
Tyz = 0.125 + 0.1*sin(2*pi*Xhr)

# density
figure(cnt.bump())
plot(Xv, num[:,4], 'ro-')
plot(Xhr, n, 'k-')
axis('tight')
xlabel('X')
ylabel('Number Density')
title('Number Density')
minorticks_on()
grid()
savefig('s7-2x3v-num.png', bbox='tight')

# momentum-x
d = gkedata.GkeData("s7-2x3v-mom_momDensity.h5")
dg1Mom = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Vc, mom = dg1Mom.project(0)

figure(cnt.bump())
plot(Xv, mom[:,4], 'ro-')
plot(Xhr, n*ux, 'k-')
xlabel('X')
ylabel('Momentum Density')
title('Momentum Density in X')
axis('tight')
minorticks_on()
grid()
savefig('s7-2x3v-momx.png', bbox='tight')

# momentum-y
Xc, Vc, mom = dg1Mom.project(1)

figure(cnt.bump())
plot(Xv, mom[:,4], 'ro-')
plot(Xhr, n*uy, 'k-')
axis('tight')
xlabel('X')
ylabel('Momentum Density')
title('Momentum Density in Y')
minorticks_on()
grid()
savefig('s7-2x3v-momy.png', bbox='tight')

# momentum-z
Xc, Vc, mom = dg1Mom.project(2)

figure(cnt.bump())
plot(Xv, mom[:,4], 'ro-')
plot(Xhr, n*uz, 'k-')
axis('tight')
xlabel('X')
ylabel('Momentum Density')
title('Momentum Density in Z')
minorticks_on()
grid()
savefig('s7-2x3v-momz.png', bbox='tight')

# total Pxx
d = gkedata.GkeData("s7-2x3v-mom_pressureTensor.h5")
dg1Pr = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Vc, pr = dg1Pr.project(0)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Txx + n*ux*ux, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{xx}$')
title(r'$P_{xx}$')
minorticks_on()
grid()
savefig('s7-2x3v-pxx.png', bbox='tight')

# total Pxy
Xc, Vc, pr = dg1Pr.project(1)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Txy + n*ux*uy, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{xy}$')
title(r'$P_{xy}$')
minorticks_on()
grid()
savefig('s7-2x3v-pxy.png', bbox='tight')

# total Pxz
Xc, Vc, pr = dg1Pr.project(2)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Txz + n*ux*uz, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{xz}$')
title(r'$P_{xz}$')
minorticks_on()
grid()
savefig('s7-2x3v-pxz.png', bbox='tight')

# total Pyy
Xc, Vc, pr = dg1Pr.project(3)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Tyy + n*uy*uy, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{yy}$')
title(r'$P_{yy}$')
minorticks_on()
grid()
savefig('s7-2x3v-pyy.png', bbox='tight')

# total Pyz
Xc, Vc, pr = dg1Pr.project(4)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Tyz + n*uy*uz, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{yz}$')
title(r'$P_{yz}$')
minorticks_on()
grid()
savefig('s7-2x3v-pyz.png', bbox='tight')

# total Pzz
Xc, Vc, pr = dg1Pr.project(5)

figure(cnt.bump())
plot(Xv, pr[:,4], 'ro-')
plot(Xhr, n*Tzz + n*uz*uz, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{zz}$')
title(r'$P_{zz}$')
minorticks_on()
grid()
savefig('s7-2x3v-pzz.png', bbox='tight')

# ptcl energy
d = gkedata.GkeData("s7-2x3v-mom_ptclEnergy.h5")
dg1Eg = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Vc, Eg = dg1Eg.project(0)

Er = 0.5*(n*(Txx+Tyy+Tzz) + n*(ux*ux+uy*uy+uz*uz))

figure(cnt.bump())
plot(Xv, Eg[:,4], 'ro-')
plot(Xhr, Er, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$Energy$')
title(r'Particle Energy')
minorticks_on()
grid()
savefig('s7-2x3v-Eg.png', bbox='tight')

show()
