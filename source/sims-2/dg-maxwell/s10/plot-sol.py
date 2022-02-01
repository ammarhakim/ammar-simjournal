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


def getRaw(q, component, numEqns, nNodes):
    rawData = numpy.zeros((q.shape[0], q.shape[1], nNodes), numpy.float)
    for n in range(nNodes):
        rawData[:,:,n] = q[:,:,component+n*numEqns]
    return rawData

def calcAverage(fld):
    wt = pylab.array([
            0.0277778,
            0.138889,
            0.138889,
            0.0277778,
            0.138889,
            0.694444,
            0.694444,
            0.138889,
            0.138889,
            0.694444,
            0.694444,
            0.138889,
            0.0277778,
            0.138889,
            0.138889,
            0.0277778,
            ])
    return (fld[:,:,0:16]*wt).sum(axis=-1)

fh = tables.openFile("s10-dg-maxwell_q_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# get final solution
q1 = getRaw(fh.root.StructGridField, 2, 8, 16)

# get intial solution
fh = tables.openFile("s10-dg-maxwell_q_0.h5")
q0 = getRaw(fh.root.StructGridField, 2, 8, 16)

# compute error
q0avg = calcAverage(q0)
q1avg = calcAverage(q1)

errAvg = 0.25*dx*dy*numpy.abs(q1avg-q0avg).sum()
print dx, errAvg

pylab.show()
