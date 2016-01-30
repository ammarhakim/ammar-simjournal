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

def calcCenters(VcEdges):
    Vx = numpy.zeros((VcEdges.shape[0]-1,), float)
    for i in range(VcEdges.shape[0]-1):
        Vx[i] = 0.5*(VcEdges[i]+VcEdges[i+1])
    return Vx

def calcExact(n, nu, E, V):
    u = nu/n
    vt2 = E-n*u**2
    return n/sqrt(2*pi*vt2)*exp(-(V-u)**2/(2*vt2))

d = gkedata.GkeData("s1-bgk-boltz_distf_0.h5")
dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
Xc, Vc, fv_0 = dg1.project(0)

d = gkedata.GkeData("s1-bgk-boltz_distf_1.h5")
dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
Xc, Vc, fv_1 = dg1.project(0)

# number density
fh = tables.openFile("s1-bgk-boltz_numDensity_0.h5")
n = fh.root.StructGridField[0,0]
# momentum density
fh = tables.openFile("s1-bgk-boltz_momentum_0.h5")
nu = fh.root.StructGridField[0,0]
# energy density
fh = tables.openFile("s1-bgk-boltz_ptclEnergy_0.h5")
E = fh.root.StructGridField[0,0]

Vx = calcCenters(Vc[:,0]) # cell center coordinates
Vex = linspace(Vc[0,0], Vc[-1,0], 32)
fex = calcExact(n, nu, E, Vex)

# plot initial and final solution
plot(Vx, fv_0[2,:], 'r-')
plot(Vx, fv_1[2,:], 'k-')
plot(Vex, fex, 'bo')

show()


