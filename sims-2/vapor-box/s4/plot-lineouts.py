from pylab import *
import tables
import euler

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

gasGamma = 5.0/3.0
amu = 1.66053892e-27 # Kg
mLi = 6.941*amu # Kg
kb = 1.3806488e-23 # J/K
Tinit = 800+273.14 # K
cs0 = sqrt(kb*Tinit/mLi)
tEnd = 5*2.0/cs0

def pressure(q):
    return euler.fluidEx.getP(q)

def mach(q):
    return euler.fluidEx.getMach(q)

def getMeshGrid(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

    return X, Y

fh = tables.openFile("s4-four-box-chain_q_10.h5")
q = fh.root.StructGridField
X, Y = getMeshGrid(fh.root.StructGrid)
nx, ny = q.shape[0], q.shape[1]

numDensity = q[:,:,0]/mLi
figure(1)
plot(X, numDensity[:,ny/2])
title('Number density [#/m^3]')
xlabel('X')
ylabel('Number Density')
axis('tight')
savefig('s4-four-box-chain-numDensity.png')

figure(2)
temp = pressure(q)/(numDensity*kb)
plot(X, temp[:,ny/2]-273.15)
title('Temperature [C]')
xlabel('X')
ylabel('Temperature')
axis('tight')
savefig('s4-four-box-chain-temperature.png')

figure(3)
machN = mach(q)
plot(X, machN[:,ny/2])
title('Mach Number')
xlabel('X')
ylabel('Mach Number')
axis('tight')
savefig('s4-four-box-chain-mach.png')

figure(4)
pr = pressure(q)
plot(X, pr[:,ny/2])
title('Pressure [Pa]')
xlabel('X')
ylabel('Pressure')
axis('tight')
savefig('s4-four-box-chain-press.png')

figure(5)
numDensity = q[:,:,0]/mLi
semilogy(X, numDensity[:,ny/2])
title('Number density [#/m^3]')
xlabel('X')
ylabel('Number Density')
axis('tight')
print(numDensity[-1,ny/2])
savefig('s4-four-box-chain-ln-numDensity.png')
