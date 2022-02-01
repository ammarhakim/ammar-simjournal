from pylab import *
import tables

import pylab
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

# function to compute lineout along ray with angle theta
def calcLineOut(nv, dx, dy, rho, theta):
    vals = zeros((nv, 1), float)
    cost = cos(theta)
    sint = sin(theta)
    dv = 1.0/nv
    vp = linspace(0.5*dv, 1.0-0.5*dv, nv)
    xp = cost*vp
    yp = sint*vp
    for i in range(nv):
        ix = (int) (xp[i]/dx)
        iy = (int) (yp[i]/dy)
        vals[i] = rho[ix,iy]
    return vals

gasGamma = 5.0/3.0

fh = tables.openFile("s400-euler-noh-ds-2d_q_5.h5")
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

figure(1)
gs = gridspec.GridSpec(1, 2, height_ratios=[1, 1])

subplot(gs[0])
# plot it
pcolormesh(XXe, YYe, transpose(rho))

rhoMin = 2.5
rhoMax = 4.0
step = 0.25
nlevels = (int) (rhoMax-rhoMin)/step
cList = linspace(rhoMin, rhoMax, nlevels)
contour(XX, YY, transpose(rho), cList, colors='k')

rhoMin = 14.0
rhoMax = 17.0
step = 0.2
nlevels = (int) (rhoMax-rhoMin)/step
cList = linspace(rhoMin, rhoMax, nlevels)
contour(XX, YY, transpose(rho), cList, colors='k')
axis('image')

# scatter plot of radial slices of density
#figure(2)
subplot(gs[1])
nv = 100
Xr = linspace(0, 1, 100)
tRange = linspace(0, pi/4, 20)
for i in range(tRange.shape[0]):
    v = calcLineOut(nv, dx, dy, rho, tRange[i])
    plot(Xr, v, '-b')

# exact solution
XrHr = linspace(0, 1, 500)
rhoEx = 0.0*XrHr
T = 2.0
for i in range(XrHr.shape[0]):
    if XrHr[i] < T/3:
        rhoEx[i] = 16.0
    else:
        rhoEx[i] = 1 + T/XrHr[i]
plot(XrHr, rhoEx, 'r-', linewidth=2)
xlabel('Radius')
ylabel('Density')

savefig('s400-noh-rho.png', dpi=300, bbox_inches='tight')
show()

