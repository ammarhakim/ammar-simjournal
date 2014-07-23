import tables
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import numpy
import math

# customization for figure
#rcParams['lines.linewidth']            = 2.0
rcParams['font.size']                  = 18
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
rcParams['figure.subplot.bottom']      = 0.125
rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

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

fh = tables.openFile("s422-euler-cyl-2d_q_3.h5")
tm = fh.root.timeData._v_attrs.vsTime
q = fh.root.StructGridField
gasGamma = 1.4

rho = q[:,:,0]
u = q[:,:,1]/rho
v = q[:,:,2]/rho
Er = q[:,:,4]
pr = (gasGamma-1)*(Er - 0.5*rho*(u*u+v*v))
cs = sqrt(gasGamma*pr/rho)
mach = sqrt(u*u+v*v)/cs

fh = tables.openFile("s422-euler-cyl-2d_inOut.h5")
inOut = fh.root.StructGridField.read()

xmin = fh.root.StructGrid._v_attrs.vsLowerBounds[0]
xmax = fh.root.StructGrid._v_attrs.vsUpperBounds[0]
NX = inOut.shape[0]
dx = (xmax-xmin)/NX
X = linspace(xmin+0.5*dx, xmax-0.5*dx, NX)

ymin = fh.root.StructGrid._v_attrs.vsLowerBounds[1]
ymax = fh.root.StructGrid._v_attrs.vsUpperBounds[1]
NY = inOut.shape[1]
dy = (ymax-ymin)/NY
Y = linspace(ymin+0.5*dy, ymax-0.5*dy, NY)

# plot density lineout
figure(1)

subplot(1,2,1)
prMa = numpy.ma.masked_where(inOut[:,:,0]<0, pr)
contour(X, Y, transpose(prMa), 20, colors='k', linewidth=0.1)
im = pcolormesh(X, Y, transpose(prMa))
title('Pressure %g' % tm)
axis('image')
colorbar_adj(im)

subplot(1,2,2)

fh = tables.openFile("s422-euler-cyl-2d_q_5.h5")
q = fh.root.StructGridField
gasGamma = 1.4

rho = q[:,:,0]
u = q[:,:,1]/rho
v = q[:,:,2]/rho
Er = q[:,:,4]
pr = (gasGamma-1)*(Er - 0.5*rho*(u*u+v*v))
cs = sqrt(gasGamma*pr/rho)
mach = sqrt(u*u+v*v)/cs

prMa = numpy.ma.masked_where(inOut[:,:,0]<0, pr)
contour(X, Y, transpose(prMa), 20, colors='k', linewidth=0.1)
im = pcolormesh(X, Y, transpose(prMa))
gca().set_yticklabels([])
tm = fh.root.timeData._v_attrs.vsTime
title('Pressure %g' % tm)
axis('image')
colorbar_adj(im)

savefig('s422-pressure.png', dpi=300)

show()
