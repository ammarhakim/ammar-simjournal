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

def colorbar_adj(obj, mode=1, redraw=False, _fig_=None, _ax_=None, aspect=None):
    '''
    Add a colorbar adjacent to obj, with a matching height
    For use of aspect, see http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.set_aspect ; E.g., to fill the rectangle, try "auto"
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if mode == 1:
        _fig_ = obj.figure; _ax_ = obj.axes
    elif mode == 2: # assume obj is in the current figure/axis instance
        _fig_ = plt.gcf(); _ax_ = plt.gca()
    _divider_ = make_axes_locatable(_ax_)
    _cax_ = _divider_.append_axes("right", size="5%", pad=0.05)
    _cbar_ =  _fig_.colorbar(obj, cax=_cax_)
    if aspect != None:
        _ax_.set_aspect(aspect)
    if redraw:
        _fig_.canvas.draw()
    return _cbar_

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
    
for i in range(0,21):
    print "Working on %d ..." % i
    d = gkedata.GkeData("r1-es-resonance_distfElc_%d.h5" % i )
    dg1 = gkedgbasis.GkeDgSerendip2DPolyOrder2Basis(d)
    Xc, Yc, fve = dg1.project(0)

    X, V = getXv(Xc, Yc)
    figure(1)
    plot(V, fve[fve.shape[0]/2, :], 'r-')
    ylo, yup = gca().get_ylim()
    plot([1.0, 1.0], [ylo, yup], '--k')
    gca().set_ylim([0, 0.4])
    title('Time %g' % d.time)
    xlabel('V')
    ylabel('f(V)')
    savefig('r1-es-resonance-fve1_%05d.png' % i)
    close()

    figure(2)
    subplot(2,1,1)
    im = pcolormesh(Xc, Yc, pylab.transpose(fve))
    plot([-2*pi,2*pi], [1.0, 1.0], 'k--', linewidth=2.0)
    axis('tight')
    colorbar_adj(im)

    subplot(2,1,2)
    im = pcolormesh(Xc, Yc, pylab.transpose(fve-fve0))
    plot([-2*pi,2*pi], [1.0, 1.0], 'k--', linewidth=2.0)
    #axis('tight')
    gca().set_ylim([0.0, 2.0])
    gca().set_xlim([-2*pi, 2*pi])
    colorbar_adj(im)
    
    suptitle('Time %g' % d.time)
    savefig('r1-es-resonance-fve_%05d.png' % i)
    close()
    
    d.close()
