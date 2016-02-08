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

def calcHamilton(Xc, Vc, mass, charge, phi):
    dx = (Xc[0,-1]-Xc[0,0])/(Xc.shape[1]-1)
    dv = (Vc[-1,0]-Vc[0,0])/(Vc.shape[0]-1)
    
    X1 = linspace(Xc[0,0]+0.5*dx, Xc[0,-1]-0.5*dx, Xc.shape[1]-1)
    V1 = linspace(Vc[0,0]+0.5*dv, Vc[-1,0]-0.5*dv, Vc.shape[0]-1)
    XX, VV = meshgrid(X1, V1)
    return XX, VV, 0.5*mass*VV**2 + charge*phi
    
for i in range(0,101):
    print "Working on %d ..." % i
    d = gkedata.GkeData("m1-es-shock_distfElc_%d.h5" % i )
    dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
    Xc, Yc, fve = dg1.project(0)

    # potential
    d = gkedata.GkeData("m1-es-shock_phi_%d.h5" % i )
    dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
    Xp, phi = dg1.project(0)
    XX, VV, ham = calcHamilton(Xc, Yc, 1.0, -1.0, phi)

    subplot(2,1,1)
    pylab.pcolormesh(Xc, Yc, pylab.transpose(fve))
    pylab.contour(XX, VV, ham, linewidth=1, levels=[0], colors='k')
    pylab.axis('tight')
    #pylab.xlabel('X')
    pylab.ylabel('V')

    d = gkedata.GkeData("m1-es-shock_distfIon_%d.h5" % i )
    dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
    Xc, Yc, fvi = dg1.project(0)
    subplot(2,1,2)
    pylab.pcolormesh(Xc, Yc, pylab.transpose(fvi))
    pylab.axis('tight')
    pylab.xlabel('X')
    pylab.ylabel('V')

    pylab.suptitle('Time %g' % d.time)
    pylab.savefig('m1-es-shock_distfElc_%05d.png' % i)
    pylab.close()
    d.close()
