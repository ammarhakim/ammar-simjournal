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

d = gkedata.GkeData("s1-em-landau_distfElc_0.h5")
dg1_0 = gkedgbasis.GkeDgSerendipNorm2DPolyOrder2Basis(d)

d = gkedata.GkeData("s1-em-landau_distfElc_1.h5")
dg1_1 = gkedgbasis.GkeDgSerendipNorm2DPolyOrder2Basis(d)

d = gkedata.GkeData("s1-em-landau_distfElc_2.h5")
dg1_2 = gkedgbasis.GkeDgSerendipNorm2DPolyOrder2Basis(d)

Xc, Yc, fve_0 = dg1_0.project(0)
Xc, Yc, fve_1 = dg1_1.project(0)
Xc, Yc, fve_2 = dg1_2.project(0)

figure(1)
subplot(2,1,1)
im = pcolormesh(Xc, Yc, pylab.transpose(fve_1-fve_0))
axis('tight')
colorbar_adj(im)

subplot(2,1,2)
im = pcolormesh(Xc, Yc, pylab.transpose(fve_2-fve_0))
axis('tight')
colorbar_adj(im)
savefig('s1-em-landau-deltaF.png')

figure(2)
im = pcolormesh(Xc, Yc, pylab.transpose(fve_2))
axis('tight')
colorbar_adj(im)
savefig('s1-em-landau-deltaF.png')
