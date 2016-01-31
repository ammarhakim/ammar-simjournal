import gkedata
import gkedgbasis
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
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

def IE(n, nu, E):
    return 0.5*(E-nu**2/n)

# density plot
d = gkedata.GkeData("../s5/s5-bgk-boltz_numDensity_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, n2 = dg1.project(0)

d = gkedata.GkeData("../s6/s6-bgk-boltz_numDensity_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, n3 = dg1.project(0)

nEul = loadtxt("../m2/m2-euler-shock-exact-density.txt")

figure(1)
plot(Xc, n2, '-r', label='Kn=1/100')
plot(Xc, n3, '-b', label='Kn=1/1000')
plot(nEul[:,0], nEul[:,1], 'k--')
xlabel('X')
ylabel('Density')
legend(loc='best')
savefig('jets-density-cmp.png', dpi=200)

# momentum plot
d = gkedata.GkeData("../s5/s5-bgk-boltz_momentum_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, nu2 = dg1.project(0)

d = gkedata.GkeData("../s6/s6-bgk-boltz_momentum_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, nu3 = dg1.project(0)

uEul = loadtxt("../m2/m2-euler-shock-exact-velocity.txt")

figure(2)
plot(Xc, nu2/n2, '-r', label='Kn=1/100')
plot(Xc, nu3/n3, '-b', label='Kn=1/1000')
plot(uEul[:,0], uEul[:,1], 'k--')
xlabel('X')
ylabel('Velocity')
legend(loc='best')
savefig('jets-velocity-cmp.png', dpi=200)

# internal energy plot
d = gkedata.GkeData("../s5/s5-bgk-boltz_ptclEnergy_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, E2 = dg1.project(0)

d = gkedata.GkeData("../s6/s6-bgk-boltz_ptclEnergy_5.h5")
dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, E3 = dg1.project(0)

pEul = loadtxt("../m2/m2-euler-shock-exact-pressure.txt")

figure(3)
plot(Xc, IE(n2, nu2, E2), '-r', label='Kn=1/100')
plot(Xc, IE(n3, nu3, E3), '-b', label='Kn=1/1000')
plot(pEul[:,0], pEul[:,1]/(3-1), 'k--')
xlabel('X')
ylabel('Particle Energy')
legend(loc='best')
savefig('jets-ptclInternalEnergy-cmp.png', dpi=200)

figure(4)
plot(Xc, 0.5*E2, '-r', label='Kn=1/100')
plot(Xc, 0.5*E3, '-b', label='Kn=1/1000')
plot(pEul[:,0], 0.5*nEul[:,1]*uEul[:,1]**2+pEul[:,1]/(3-1), 'k--')
xlabel('X')
ylabel('Particle Energy')
legend(loc='best')
savefig('jets-ptclEnergy-cmp.png', dpi=200)

show()


