from pylab import *
import tables

import pylab
from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
#rcParams['lines.linewidth']            = 2
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

gasGamma = 1.4

fh = tables.openFile("s398-euler-reim-ds-2d_q_5.h5")
grid = fh.root.StructGrid
nx, ny = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1]
dx = 1.0/nx
dy = 1.0/ny
Xe = linspace(0.0, 1.0, nx+1)
Ye = linspace(0.0, 1.0, ny+1)
XXe, YYe = meshgrid(Xe, Ye)

X = linspace(0.5*dx, 1.0-0.5*dx, nx)
Y = linspace(0.5*dx, 1.0-0.5*dx, ny)
XX, YY = meshgrid(X, Y)

q = fh.root.StructGridField
rho = q[:,:,0]
u = q[:,:,1]/rho
v = q[:,:,2]/rho
pr = (q[:,:,4] - 0.5*rho*(u**2+v**2))*(gasGamma-1)

# plot it
pcolormesh(XXe, YYe, transpose(pr))
colorbar()

rhoMin = 0.43
rhoMax = 0.99
step = 0.02
nlevels = (int) (rhoMax-rhoMin)/step
cList = linspace(rhoMin, rhoMax, 29)
contour(XX, YY, transpose(rho), cList, colors='k')

savefig('s398-pr-dens-flow.png', dpi=300, bbox_inches='tight')

show()
