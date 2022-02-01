from pylab import *

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

domsz = array([5.0, 10.0, 15.0, 25.0])
rates = array([0.9583, 0.8196, 0.7546, 0.6766])

# compute coefficients
f1 = (rates[1]-rates[0])/(1/sqrt(domsz[1])-1/sqrt(domsz[0]))
f0 = rates[0] - f1/sqrt(domsz[0])

# compute smooth curve
dhr = linspace(5.0, 100.0, 50)
rsqrt = f0 + f1/sqrt(dhr)

# now make plots
plot(domsz, rates, 'ro', linewidth=2, label='10M')
plot(dhr, rsqrt, 'k-', label=r'$1/sqrt(di/L)$')
title('Scaling of 10-moment peak reconnection rates')
xlabel('Island size')
ylabel('Peak Rate')
legend()
savefig('peak-rates-10M.png')
show()

