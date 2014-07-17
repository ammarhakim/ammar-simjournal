import tables
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import numpy
import math

# customization for figure
rcParams['lines.linewidth']            = 2.0
rcParams['font.size']                  = 18
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
rcParams['figure.subplot.bottom']      = 0.125
rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

def getVars(q):
    rho = q[:,:,0]
    u = q[:,:,1]/rho
    v = q[:,:,2]/rho
    Er = q[:,:,4]
    pr = (gasGamma-1)*(Er - 0.5*rho*(u*u+v*v))
    cs = sqrt(gasGamma*pr/rho)
    mach = sqrt(u*u+v*v)/cs
    return rho, pr, mach
    

fh = tables.openFile("s417-euler-wedge-2d_q_10.h5")
q = fh.root.StructGridField
gasGamma = 1.4

rho_s417, pr_s417, mach_s417 = getVars(q)

xmin = fh.root.StructGrid._v_attrs.vsLowerBounds[0]
xmax = fh.root.StructGrid._v_attrs.vsUpperBounds[0]
NX = fh.root.StructGrid._v_attrs.vsNumCells[0]
dx = (xmax-xmin)/NX
X_s417 = linspace(xmin+0.5*dx, xmax-0.5*dx, NX)

ymin = fh.root.StructGrid._v_attrs.vsLowerBounds[1]
ymax = fh.root.StructGrid._v_attrs.vsUpperBounds[1]
NY = fh.root.StructGrid._v_attrs.vsNumCells[1]
dy = (ymax-ymin)/NY
Y_s417 = linspace(ymin+0.5*dy, ymax-0.5*dy, NY)

fh = tables.openFile("s417-euler-wedge-2d_inOut.h5")
inOut_s417 = fh.root.StructGridField.read()

##

fh = tables.openFile("../s416/s416-euler-wedge-2d_q_10.h5")
q = fh.root.StructGridField
gasGamma = 1.4

rho_s416, pr_s416, mach_s416 = getVars(q)

fh = tables.openFile("../s416/s416-euler-wedge-2d_inOut.h5")
inOut_s416 = fh.root.StructGridField.read()

xmin = fh.root.StructGrid._v_attrs.vsLowerBounds[0]
xmax = fh.root.StructGrid._v_attrs.vsUpperBounds[0]
NX = fh.root.StructGrid._v_attrs.vsNumCells[0]
dx = (xmax-xmin)/NX
X_s416 = linspace(xmin+0.5*dx, xmax-0.5*dx, NX)

ymin = fh.root.StructGrid._v_attrs.vsLowerBounds[1]
ymax = fh.root.StructGrid._v_attrs.vsUpperBounds[1]
NY = fh.root.StructGrid._v_attrs.vsNumCells[1]
dy = (ymax-ymin)/NY
Y_s416 = linspace(ymin+0.5*dy, ymax-0.5*dy, NY)

# 2D plots of density pressure and Mach number
figure(1)

# density
subplot(3,2,1)
pcolor(X_s416, Y_s416, transpose(numpy.ma.masked_where(inOut_s416[:,:,0]<0, rho_s416)))
title('Density')
ylabel('Y')
gca().get_xaxis().set_ticklabels([])
axis('image')

subplot(3,2,2)
pcolor(X_s417, Y_s417, transpose(numpy.ma.masked_where(inOut_s417[:,:,0]<0, rho_s417)))
title('Density')
gca().get_xaxis().set_ticklabels([])
gca().get_yaxis().set_ticklabels([])
axis('image')

# pressure
subplot(3,2,3)
pcolor(X_s416, Y_s416, transpose(numpy.ma.masked_where(inOut_s416[:,:,0]<0, pr_s416)))
title('Pressure')
ylabel('Y')
gca().get_xaxis().set_ticklabels([])
axis('image')

subplot(3,2,4)
pcolor(X_s417, Y_s417, transpose(numpy.ma.masked_where(inOut_s417[:,:,0]<0, pr_s417)))
title('Pressure')
gca().get_xaxis().set_ticklabels([])
gca().get_yaxis().set_ticklabels([])
axis('image')

# mach
subplot(3,2,5)
pcolor(X_s416, Y_s416, transpose(numpy.ma.masked_where(inOut_s416[:,:,0]<0, mach_s416)))
title('Mach Number')
ylabel('Y')
xlabel('X')
axis('image')

subplot(3,2,6)
pcolor(X_s417, Y_s417, transpose(numpy.ma.masked_where(inOut_s417[:,:,0]<0, mach_s417)))
title('Mach Number')
xlabel('X')
gca().get_yaxis().set_ticklabels([])
axis('image')

savefig('s416-417-wedge-cmp.png', dpi=300)
show()
