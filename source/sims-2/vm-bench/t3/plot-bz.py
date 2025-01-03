from pylab import *
import tables

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

nFrame = 75
for i in range(0, nFrame+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("t3-5m-wiebel_q_%d.h5" % i)
    q = fh.root.StructGridField
    Bz = q[:,:,15]
    gd = fh.root.StructGrid._v_attrs
    Y = linspace(gd.vsLowerBounds[1], gd.vsUpperBounds[1], gd.vsNumCells[1])

    figure(1)
    plot(Y, Bz[0,:], 'k-')
    tm = float(fh.root.timeData._v_attrs.vsTime)
    title('T=%g' % float(fh.root.timeData._v_attrs.vsTime))
    axis('tight')    
    gca().set_ylim([-0.15,0.15])
    savefig('t3-5m-wiebel_Bz_%05d.png' % i)
    close()

    fh.close()
