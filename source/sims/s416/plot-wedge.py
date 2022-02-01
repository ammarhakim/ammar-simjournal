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
rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
rcParams['figure.subplot.bottom']      = 0.125
rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

fh = tables.openFile("s416-euler-wedge-2d_q_10.h5")
q = fh.root.StructGridField
gasGamma = 1.4

rho = q[:,:,0]
u = q[:,:,1]/rho
v = q[:,:,2]/rho
Er = q[:,:,4]
pr = (gasGamma-1)*(Er - 0.5*rho*(u*u+v*v))
cs = sqrt(gasGamma*pr/rho)
mach = sqrt(u*u+v*v)/cs

fh = tables.openFile("s416-euler-wedge-2d_inOut.h5")
inOut = fh.root.StructGridField.read()

xmin = fh.root.StructGrid._v_attrs.vsLowerBounds[0]
xmax = fh.root.StructGrid._v_attrs.vsUpperBounds[0]
NX = inOut.shape[0]
dx = (xmax-xmin)/NX

ymin = fh.root.StructGrid._v_attrs.vsLowerBounds[1]
ymax = fh.root.StructGrid._v_attrs.vsUpperBounds[1]
NY = inOut.shape[1]
dy = (ymax-ymin)/NY
Y = linspace(ymin+0.5*dy, ymax-0.5*dy, NY)

idx = int(0.9/dx)

# plot density lineout
figure(1)
subplot(1,2,1)
rhoMa = numpy.ma.masked_where(inOut[:,:,0]<0, rho)
plot(rhoMa[idx,:], Y, '-m')
plot([3.7126, 3.7126], [0.2, 0.55], 'k--', linewidth=1.0)
gca().get_xaxis().set_ticks([1.0, 2.0, 3.0, 4.0])
gca().set_ylim([0.2,0.55])
title('Density at x=0.9')
xlabel('Density')
ylabel('Y')

subplot(1,2,2)
prMa = numpy.ma.masked_where(inOut[:,:,0]<0, pr)
plot(prMa[idx,:], Y, '-m')
plot([9.3, 9.3], [0.2, 0.55], 'k--', linewidth=1.0)
gca().get_yaxis().set_ticks([])
gca().set_ylim([0.2,0.55])
title('Pressure at x=0.9')
xlabel('Pressure')

savefig('s416-density-pressure-lineout.png')

figure(3)
machMa = numpy.ma.masked_where(inOut[:,:,0]<0, mach)
plot(machMa[idx,:], Y, '-m')
title('Vertical lineout of Mach number')
xlabel('Mach Number')
ylabel('Y')

show()
