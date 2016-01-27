import gkedata
import gkedgbasis

import pylab
import tables
import math
import numpy
pylab.rc('text', usetex=True)

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

NT = 100
Ey_tx = numpy.zeros( (NT+1, 64*3), numpy.float )
for i in range(NT+1):
    dg1 = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(gkedata.GkeData("t1-plasma-beach_em_%d.h5" % i))
    X, Ey = dg1.project(1)
    Ey_tx[i,:] = Ey

dx = 1/400.0
T = pylab.linspace(0, 5e-9, NT+1)
TT, XX = pylab.meshgrid(T, X)

# compute cutoff location
dx100 = 1/100.
deltaT = dx100/2.99792458e8
driveOmega = 3.14159265358979323846264338328/10/deltaT
xcutoff = 1-math.pow(driveOmega/25*deltaT, 1/5.0)

pylab.pcolormesh(TT, XX, Ey_tx.transpose())
pylab.plot([0, 5e-9], [xcutoff, xcutoff], 'k--', linewidth=2)
pylab.xlabel('Time [s]')
pylab.ylabel('X [m]')
pylab.title(r'$E_y(t,x)$')
pylab.savefig('t1-plasmabeach_Ey.png')
pylab.show()
