import pylab
import tables
import math
import numpy
import pylab
import numpy
import gkedata
import gkedgbasis
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

def plotFrame(f):
    gkd = gkedata.GkeData("s2-dg-euler-rt_q_%d.h5" % f)
    dgd = gkedgbasis.GkeDgSerendip2DPolyOrder1Basis(gkd)
    Xc, Yc, rho = dgd.project(0)
    pylab.figure(1)
    pylab.pcolormesh(Xc, Yc, pylab.transpose(rho))
    pylab.axis('image')
    pylab.gca().set_xticks([])
    pylab.gca().set_yticks([])
    pylab.axis('image')
    pylab.clim(35000, 65000)
    pylab.title ('t = %g' % gkd.time)
    pylab.savefig('s2-dg-euler-rt_rho_%05d.png' % f, bbox_inches='tight')
    pylab.close()

for i in range(0,26):
    print ("Working on frame %d ..." % i)
    plotFrame(i)
