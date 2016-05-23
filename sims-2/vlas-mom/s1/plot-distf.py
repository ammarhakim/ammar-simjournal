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

def getXv(Xc, Vc):
    dx = (Xc[0,-1]-Xc[0,0])/(Xc.shape[1]-1)
    dv = (Vc[-1,0]-Vc[0,0])/(Vc.shape[0]-1)

    X1 = linspace(Xc[0,0]+0.5*dx, Xc[0,-1]-0.5*dx, Xc.shape[1]-1)
    V1 = linspace(Vc[0,0]+0.5*dv, Vc[-1,0]-0.5*dv, Vc.shape[0]-1)

    return X1, V1

# density
d = gkedata.GkeData("s1-1x1v-mom_distf.h5")
dg1Num = gkedgbasis.GkeDgSerendipNorm2DPolyOrder1Basis(d)
Xc, Yc, fve = dg1Num.project(0)

figure(1)
pcolormesh(Xc, Yc, transpose(fve))

show()



