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

def evalSum(coeff, fields):
    res = 0.0*fields[0]
    for i in range(len(coeff)):
        res = res + coeff[i]*fields[i]
    return res

def projectOnFinerGrid_f(Xc, Yc, q):
    dx = Xc[1]-Xc[0]
    dy = Yc[1]-Yc[0]
    nx = Xc.shape[0]
    ny = Yc.shape[0]

    # mesh coordinates
    Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 2*nx+1) # one more
    Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 2*ny+1) # one more
    XXn, YYn = pylab.meshgrid(Xn, Yn)

    # data
    qn = pylab.zeros((2*Xc.shape[0], 2*Yc.shape[0]), float)

    v1 = q[:,:,0]
    v2 = q[:,:,1]
    v3 = q[:,:,2]
    v4 = q[:,:,3]

    vList = [v1,v2,v3,v4]

    # node 1
    c1 = [0.5625,0.1875,0.0625,0.1875]
    qn[0:2*nx:2, 0:2*ny:2] = evalSum(c1, vList)

    # node 2
    c2 = [0.1875,0.5625,0.1875,0.0625]
    qn[1:2*nx:2, 0:2*ny:2] = evalSum(c2, vList)

    # node 3
    c3 = [0.1875,0.0625,0.1875,0.5625]
    qn[0:2*nx:2, 1:2*ny:2] = evalSum(c3, vList)

    # node 4
    c4 = [0.0625,0.1875,0.5625,0.1875]
    qn[1:2*nx:2, 1:2*ny:2] = evalSum(c4, vList)
   
    return XXn, YYn, qn

fh = tables.openFile("s3-dg-maxwell_q_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# get final solution
q = getRaw(fh.root.StructGridField, 2, 8, 4)
Xn, Yn, qn_1 = projectOnFinerGrid_f(Xc, Yc, q)

# get intial solution
fh = tables.openFile("s3-dg-maxwell_q_0.h5")
q = getRaw(fh.root.StructGridField, 2, 8, 4)
Xn, Yn, qn_0 = projectOnFinerGrid_f(Xc, Yc, q)

nx, ny = Xn.shape[0], Yn.shape[0]

# make plot
pylab.figure(1)

pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1))
pylab.axis('tight')
pylab.savefig('s3-dg-maxwell-Ez.png')

# compute error
err = numpy.abs(qn_1-qn_0).sum()/(nx*ny)
print math.sqrt(dx*dy), err

pylab.show()
