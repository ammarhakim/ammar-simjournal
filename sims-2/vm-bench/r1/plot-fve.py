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

# initial conditions
d = gkedata.GkeData("r1-es-resonance_distfElc_0.h5")
dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
Xc, Yc, fve0 = dg1.project(0)
    
for i in range(0,11):
    print "Working on %d ..." % i
    d = gkedata.GkeData("r1-es-resonance_distfElc_%d.h5" % i )
    dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
    Xc, Yc, fve = dg1.project(0)

    X, V = getXv(Xc, Yc)
    figure(1)
    plot(V, fve[fve.shape[0]/2, :], 'r-')
    title('Time %g' % d.time)
    xlabel('V')
    ylabel('f(V)')
    savefig('r1-es-resonance-fve1_%05d.png' % i)
    close()

    figure(2)
    pcolormesh(Xc, Yc, pylab.transpose(fve-fve0))
    title('Time %g' % d.time)
    colorbar()
    axis('tight')
    savefig('r1-es-resonance-fve_%05d.png' % i)
    close()    
    
    d.close()
