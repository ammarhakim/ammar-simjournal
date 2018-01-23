import postgkyl
import postgkyl
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

def getXv(Xc, Vc):
    dx = (Xc[0,-1]-Xc[0,0])/(Xc.shape[1]-1)
    dv = (Vc[-1,0]-Vc[0,0])/(Vc.shape[0]-1)

    X1 = linspace(Xc[0,0]+0.5*dx, Xc[0,-1]-0.5*dx, Xc.shape[1]-1)
    V1 = linspace(Vc[0,0]+0.5*dv, Vc[-1,0]-0.5*dv, Vc.shape[0]-1)

    return X1, V1

# density
d = postgkyl.GData("x1-1x1v-max-mom_numDensity.bp")
dg1Num = postgkyl.GInterpModal(d, 1, "mo")
Xc, num = dg1Num.interpolate(0)

Xhr = linspace(Xc[0][0], Xc[0][-1], 200) # for plotting

n = sin(2*pi*Xhr)
ux = 0.1*cos(2*pi*Xhr)
Txx = 0.75 + 0.25*cos(2*pi*Xhr)

# density
figure(1)
plot(Xc[0], num, 'ro-')
plot(Xhr, n, 'k-')
axis('tight')
xlabel('X')
ylabel('Number Density')
title('Number Density')
minorticks_on()
grid()
savefig('x1-1x1v-max-num.png', bbox='tight')

# momentum
d = postgkyl.GData("x1-1x1v-max-mom_momentum.bp")
dg1Mom = postgkyl.GInterpModal(d, 1, "mo")
Xc, mom = dg1Mom.interpolate(0)

figure(2)
plot(Xc[0], mom, 'ro-')
plot(Xhr, n*ux, 'k-')
axis('tight')
xlabel('X')
ylabel('Momentum Density')
title('Momentum Density')
minorticks_on()
grid()
savefig('x1-1x1v-max-mom.png', bbox='tight')

# total Pxx
d = postgkyl.GData("x1-1x1v-max-mom_pressureTensor.bp")
dg1Pr = postgkyl.GInterpModal(d, 1, "mo")
Xc, pr = dg1Pr.interpolate(0)

figure(3)
plot(Xc[0], pr, 'ro-')
plot(Xhr, n*Txx + n*ux*ux, 'k-')
axis('tight')
xlabel('X')
ylabel(r'$P_{xx}$')
title(r'$P_{xx}$')
minorticks_on()
grid()
savefig('x1-1x1v-max-pr.png', bbox='tight')

# ptcl energy
d = postgkyl.GData("x1-1x1v-max-mom_ptclEnergy.bp")
dg1Eg = postgkyl.GInterpModal(d, 1, "mo")
Xc, Eg = dg1Eg.interpolate(0)

figure(4)
plot(Xc[0], Eg, 'ro-')
plot(Xhr, (n*Txx + n*ux*ux), 'k-')
axis('tight')
xlabel('X')
ylabel(r'$Energy$')
title(r'Particle Energy')
minorticks_on()
grid()
savefig('x1-1x1v-max-Eg.png', bbox='tight')

show()
