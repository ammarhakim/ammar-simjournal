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

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p', '--plot', action = 'store',
                  dest = 'fileName',
                  help = 'Hdf5 file to plot')

(options, args) = parser.parse_args()
fileName = options.fileName

gasGamma = 1.4

fh = tables.open_file(fileName)
grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx = lx/nx
dy = ly/ny
dz = lz/nz
Xe = linspace(0.0, lx, nx+1)
Ye = linspace(0.0, ly, ny+1)
Ze = linspace(0.0, lz, nz+1)

Xc = linspace(0.5*dx, lx-0.5*dx, nx)
Yc = linspace(0.5*dy, ly-0.5*dy, ny)
Zc = linspace(0.5*dz, ly-0.5*dz, ny)

q = fh.root.StructGridField
rho = q[:,:,:,0]
u = q[:,:,:,1]/rho
v = q[:,:,:,2]/rho
pr = (q[:,:,:,4] - 0.5*rho*(u**2+v**2))*(gasGamma-1)

# plot it
figure(1)
pcolormesh(Xe, Ze, transpose(rho[:,0,:]))
colorbar()
title('Y=0')
axis('image')

figure(2)
pcolormesh(Xe, Ye, transpose(rho[:,:,-1]))
colorbar()
title('Z=1')
axis('image')

figure(3)
pcolormesh(Xe, Ye, transpose(rho[:,:,0]))
colorbar()
title('Z=0')
axis('image')

show()
