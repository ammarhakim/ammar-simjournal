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

d = gkedata.GkeData("s1-bgk-boltz_distf_0.h5")
dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
Xc, Vc, fv_0 = dg1.project(0)

d = gkedata.GkeData("s1-bgk-boltz_distf_1.h5")
dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
Xc, Vc, fv_1 = dg1.project(0)

Vext = linspace(-6, 6, 96)

# plot initial and final solution
plot(Vext, fv_0[2,:], 'r-')
plot(Vext, fv_1[2,:], 'k-')

show()


