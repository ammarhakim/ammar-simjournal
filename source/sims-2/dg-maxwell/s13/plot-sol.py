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
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
#rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

figure(1)

subplot(1,2,1)
dg1 = gkedgbasis.GkeDgLobatto2DPolyOrder3Basis(gkedata.GkeData("s13-dg-maxwell_q_1.h5"))
Xc, Yc, Ez1 = dg1.project(2)
pcolormesh(Xc, Yc, transpose(Ez1))
title('t=1.5')
axis('image')

subplot(1,2,2)
dg2 = gkedgbasis.GkeDgLobatto2DPolyOrder3Basis(gkedata.GkeData("s13-dg-maxwell_q_2.h5"))
Xc, Yc, Ez2 = dg2.project(2)
pcolormesh(Xc, Yc, transpose(Ez2))
title('t=3.0')
axis('image')

pylab.savefig('s13-ez.png', bbox_inches='tight')

show()

