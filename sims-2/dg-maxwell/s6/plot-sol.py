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
    Xn = pylab.linspace(Xc[0]-0.5*dx, Xc[-1]+0.5*dx, 3*nx+1) # one more
    Yn = pylab.linspace(Yc[0]-0.5*dy, Yc[-1]+0.5*dy, 3*ny+1) # one more
    XXn, YYn = pylab.meshgrid(Xn, Yn)

    # data
    qn = pylab.zeros((3*Xc.shape[0], 3*Yc.shape[0]), float)

    v1 = q[:,:,0]
    v2 = q[:,:,1]
    v3 = q[:,:,2]
    v4 = q[:,:,3]
    v5 = q[:,:,4]
    v6 = q[:,:,5]
    v7 = q[:,:,6]
    v8 = q[:,:,7]

    vList = [v1,v2,v3,v4,v5,v6,v7,v8]

    # node 1
    c1 = [.2314814814814815,-.1388888888888889,-.06481481481481481,-.1388888888888889,0.462962962962963,.09259259259259259,.09259259259259259,0.462962962962963]
    qn[0:3*nx:3, 0:3*ny:3] = evalSum(c1, vList)

    # node 2
    c2 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.8333333333333334,.2777777777777778,.1666666666666667,.2777777777777778]
    qn[1:3*nx:3, 0:3*ny:3] = evalSum(c2, vList)

    # node 3
    c3 = [-.1388888888888889,.2314814814814815,-.1388888888888889,-.06481481481481481,0.462962962962963,0.462962962962963,.09259259259259259,.09259259259259259]
    qn[2:3*nx:3, 0:3*ny:3] = evalSum(c3, vList)

    # node 4
    c4 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.2777777777777778,.1666666666666667,.2777777777777778,.8333333333333334]
    qn[0:3*nx:3, 1:3*ny:3] = evalSum(c4, vList)

    # node 5
    c5 = [-0.25,-0.25,-0.25,-0.25,0.5,0.5,0.5,0.5]
    qn[1:3*nx:3, 1:3*ny:3] = evalSum(c5, vList)

    # node 6
    c6 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.2777777777777778,.8333333333333334,.2777777777777778,.1666666666666667]
    qn[2:3*nx:3, 1:3*ny:3] = evalSum(c6, vList)

    # node 7
    c7 = [-.1388888888888889,-.06481481481481481,-.1388888888888889,.2314814814814815,.09259259259259259,.09259259259259259,0.462962962962963,0.462962962962963]
    qn[0:3*nx:3, 2:3*ny:3] = evalSum(c7, vList)

    # node 8
    c8 = [-.1388888888888889,-.1388888888888889,-.1388888888888889,-.1388888888888889,.1666666666666667,.2777777777777778,.8333333333333334,.2777777777777778]
    qn[1:3*nx:3, 2:3*ny:3] = evalSum(c8, vList)

    # node 9
    c9 = [-.06481481481481481,-.1388888888888889,.2314814814814815,-.1388888888888889,.09259259259259259,0.462962962962963,0.462962962962963,.09259259259259259]
    qn[2:3*nx:3, 2:3*ny:3] = evalSum(c9, vList)
   
    return XXn, YYn, qn            

fh = tables.openFile("s6-dg-maxwell_q_1.h5")
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

dx = (upper[0]-lower[0])/cells[0]
dy = (upper[1]-lower[1])/cells[1]

Xc = pylab.linspace(lower[0]+0.5*dx, upper[0]-0.5*dx, cells[0])
Yc = pylab.linspace(lower[1]+0.5*dy, upper[1]-0.5*dy, cells[1])

# get final solution
q1 = getRaw(fh.root.StructGridField, 2, 8, 8)
Xn, Yn, qn_1 = projectOnFinerGrid_f(Xc, Yc, q1)

# get intial solution
fh = tables.openFile("s6-dg-maxwell_q_0.h5")
q0 = getRaw(fh.root.StructGridField, 2, 8, 8)
Xn, Yn, qn_0 = projectOnFinerGrid_f(Xc, Yc, q0)

nx, ny = Xn.shape[0], Yn.shape[0]

# make plot
pylab.figure(1)

pylab.pcolormesh(Xn, Yn, pylab.transpose(qn_1))
pylab.axis('tight')
pylab.savefig('s6-dg-maxwell-Ez.png')

def calcAverage(fld):
    d13 = 1.0/3.0
    d43 = 4.0/3.0
    wt = pylab.array([-d13, -d13, -d13, -d13, d43, d43, d43, d43])
    return (fld[:,:,0:8]*wt).sum(axis=-1)

# compute error
q0avg = calcAverage(q0)
q1avg = calcAverage(q1)
vol = dx*dy/4.0

errAvg = vol*numpy.abs(q1avg-q0avg).sum()
print dx, errAvg

pylab.show()
