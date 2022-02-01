import numpy
import pylab
import tables
import optparse
import gkedata
import gkedgbasis

from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
#rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 14
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 1 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

gasGamma = 1.4

# read data
d = gkedata.GkeData("s4-dg-euler-pos_q_1.h5")
dgDat = gkedgbasis.GkeDgLobatto1DPolyOrder1Basis(d)
X, rho = dgDat.project(0)
X, rhou = dgDat.project(1)
X, Er = dgDat.project(4)

u = rhou/rho
pr = (gasGamma-1)*(Er - 0.5*rho*u*u)

ex_density = pylab.loadtxt("s27-euler-shock-exact-density.txt")
ex_velocity = pylab.loadtxt("s27-euler-shock-exact-velocity.txt")
ex_pressure = pylab.loadtxt("s27-euler-shock-exact-pressure.txt")
ex_ie = pylab.loadtxt("s27-euler-shock-exact-internal-energy.txt")

pylab.figure(1)
pylab.subplot(2, 2, 1)
pylab.plot(X, rho, 'k-')
pylab.plot(ex_density[:,0], ex_density[:,1], 'r-')
pylab.axis('tight')
pylab.ylabel("Density")

pylab.subplot(2, 2, 2)
pylab.plot(X, u, 'k-')
pylab.plot(ex_velocity[:,0], ex_velocity[:,1], 'r-')
pylab.axis('tight')
pylab.ylabel("Velocity")

pylab.subplot(2, 2, 3)
pylab.plot(X, pr, 'k-')
pylab.plot(ex_pressure[:,0], ex_pressure[:,1], 'r-')
pylab.axis('tight')
pylab.ylabel("Pressure")

pylab.subplot(2, 2, 4)
pylab.plot(X, pr/rho/(gasGamma-1), 'k-')
pylab.plot(ex_ie[:,0], ex_ie[:,1], 'r-')
pylab.axis('tight')
pylab.ylabel("Internal Energy")

pylab.savefig("s4-dg-euler-pos_exact_cmp.png")#, bbox_inches='tight')

pylab.show()

    
