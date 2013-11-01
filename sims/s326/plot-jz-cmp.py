from pylab import *
import tables

import pylab
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

Lx = 8*pi
Ly = 4*pi

X = linspace(-Lx/2, Lx/2, 1001)
Y = linspace(-Ly/2, Ly/2, 501)

fh = tables.openFile("../s314/s314-5m-gem_q_120.h5")
q = fh.root.StructGridField
jzeAda = -25*q[:,:,3]

fh = tables.openFile("../s326/s326-iso-gem_q_24.h5")
q = fh.root.StructGridField
jzeIso = -25*q[:,:,3]

subplot(2,1,1)
im = pcolormesh(X, Y, transpose(jzeAda))
axis('image')
colorbar_adj(im)

subplot(2,1,2)
im = pcolormesh(X, Y, transpose(jzeIso))
axis('image')
colorbar_adj(im)
suptitle('Ideal Gas (top) v/s Isothermal Fluid (bottom)')

savefig('s326-s314-jze-cmp.png')

figure(2)
dat = loadtxt('s326-iso-gem_byFlux.txt')
tIso = dat[:,0]
byFluxIso = dat[:,1]

dat = loadtxt('../s314/s314-5m-gem_byFlux.txt')
tAda = dat[:,0]
byFluxAda = dat[:,1]

plot(tIso*0.1, byFluxIso, 'r-', label='Isothermal')
plot(tAda*0.1, byFluxAda, 'k-', label='Ideal gas')
xlabel('Time')
ylabel('Reconnected Flux')
legend(loc='best')

savefig('s326-314-byFlux-cmp.png')

show()


