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

d = gkedata.GkeData("2x2v-pij_numDensity.h5")
dg1Num = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Yc, num = dg1Num.project(0)

d = gkedata.GkeData("2x2v-pij_momDensity.h5")
dg1Mom = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Yc, mom = dg1Mom.project(0)

# compute velocity
#ux = dg1Mom.q[:,:,:]/dg1Num.q[:,:,:]
#uxInterp = gkedgbasis.interpOnMesh2D(gkedginterpdat.GkeDgSerendipNorm2DPolyOrder1Basis.cMat_i2 , ux)

X1, V1 = getXv(Xc, Yc) # 1D cell-centers

# make plots
figure(1)
plot(X1, num[:,4], 'ro-')
plot(X1, sin(2*pi*X1), 'k-')
xlabel('X')
ylabel('Number Density')
title('Number Density')
savefig('2x2v-num.png', bbox='tight')

figure(2)
plot(X1, mom[:,4], 'ro-')
plot(X1, 0.1*cos(2*pi*X1)*sin(2*pi*X1), 'k-')
xlabel('X')
ylabel(r'$nu_x$')
title(r'X-Momentum')
savefig('2x2v-ux.png', bbox='tight')

show()



