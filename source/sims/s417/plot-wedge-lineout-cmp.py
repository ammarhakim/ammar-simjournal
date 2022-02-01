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

idx_s417 = int(0.9/dx)

##

ys = 0.9*math.tan(20.9*pi/180)

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

idx_s416 = int(0.9/dx)

# 1D lineouts
figure(1)

# density

subplot(1,2,1)

rho = numpy.ma.masked_where(inOut_s417[:,:,0]<0, rho_s417)
plot(rho[idx_s417,:], Y_s417, '-r')
rho = numpy.ma.masked_where(inOut_s416[:,:,0]<0, rho_s416)
plot(rho[idx_s416,:], Y_s416, '-k')
plot([3.7126], [ys], 'mo', markersize=10)
plot([3.7126], [ys], 'kx')
gca().get_xaxis().set_ticks([1,2,3,4])
axis('tight')
xlabel('Density')
ylabel('Y')

subplot(1,2,2)

pr = numpy.ma.masked_where(inOut_s417[:,:,0]<0, pr_s417)
plot(pr[idx_s417,:], Y_s417, '-r')
pr = numpy.ma.masked_where(inOut_s416[:,:,0]<0, pr_s416)
plot(pr[idx_s416,:], Y_s416, '-k')
plot([9.3], [ys], 'mo', markersize=10)
plot([9.3], [ys], 'kx')
axis('tight')
xlabel('Pressure')
gca().get_yaxis().set_ticklabels([])

savefig('s416-417-wedge-lineout-cmp.png', dpi=300)
show()
