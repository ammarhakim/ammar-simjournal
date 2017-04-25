from scipy import special
import math
from pylab import *
import numpy.ma

import pylab as plt
plt.style.use('/Users/ahakim/.config/matplotlib/stylelib/pplot.mplstyle')

from optparse import OptionParser

# set command line options
parser = OptionParser()
parser.add_option('-k', '--wnumber', action = 'store',
                  dest = 'kInp', help = 'Wave number')
parser.add_option('-r', '--wR', action = 'store',
                  dest = 'wR', help = 'Maximum real extent')
parser.add_option('-T', '--TiTe', action = 'store',
                  dest = 'TiTe', help = 'Ti/Te')
parser.add_option('-i', '--wI', action = 'store',
                  dest = 'wI', help = 'Maximum imag extent')

parser.add_option('-f', action = 'store',
                  dest = 'f', help = 'Maximum imag extent')

(options, args) = parser.parse_args()

def plasmaDisp(z):
    r"""Plasma dispersion function Z(z) expressed in terms of complex error
    function.
    """
    return 1j*sqrt(pi)*exp(-z**2)*(1+special.erf(1j*z))

def derivPlasmaDisp(z):
    r"""Derivative of the plasma dispersion function.
    """
    return -2*(1+z*plasmaDisp(z))

# Some parameters
mi = 1.0
mu = 1/1836.2
vthi = 1.0
wpi = 1.0
Ti_Te = float(options.TiTe)
gasGamma = 3.0
k = float(options.kInp)

Ti = mi*vthi**2
Te = Ti/Ti_Te

me = mi*mu

wpe = wpi/math.sqrt(mu)
vthe = sqrt(Te/me)
cse = sqrt(gasGamma)*vthi/math.sqrt(Ti_Te)/math.sqrt(mu)
csi = sqrt(3)*vthi

wwci = sqrt(wpi**2 + k**2*csi**2)
print("Plasma frequency = %g. Electron sound speed = %g" % (wpe, cse))
print("Approximate real frequency %g" % wwci)

def iawEps(w, k):
    return 1-wpi**2/(2*vthi**2*k**2)*derivPlasmaDisp(w/(math.sqrt(2)*vthi*k)) - wpe**2/(w**2-cse**2*k**2)

def iawKinEps(w, k):
    return 1-wpi**2/(2*vthi**2*k**2)*derivPlasmaDisp(w/(math.sqrt(2)*vthi*k)) - wpe**2/(2*vthe**2*k**2)*derivPlasmaDisp(w/(math.sqrt(2)*vthe*k))

def iawFluid(w, k):
    return 1-wpi**2/(w**2-csi**2*k**2) - wpe**2/(w**2-cse**2*k**2)

wRl = options.wR.split(",")
wR_l = float(wRl[0])
wR_u = float(wRl[1])

wIl = options.wI.split(",")
wI_l = float(wIl[0])
wI_u = float(wIl[1])

X = linspace(wR_l, wR_u, 400)
Y = linspace(wI_l, wI_u, 400)
XX, YY = meshgrid(X, Y)
ZZ = XX+1j*YY

figure(1)
iawDisp = iawEps(ZZ, k)
iawKinDisp = iawKinEps(ZZ, k)

# The trick here is to plot the zero contours of the real-part of the
# dispersion function, and the zero contours of the imaginary-part of
# the dispersion function. The intersection of these two curves gives
# the roots.
#
# There will be more than one root. For damped modes, we need to look
# at roots in the lower half-plane which are closest to the Imag(Z)=0
# line. For growing modes, we need to look at roots in the upper
# half-plane which are farthest to the Imag(Z)=0 line.

figure(1)
suptitle('Hybrid-Kinetic (Left) Full-Kinetic (Right) k*vti/wpi=%g, Ti/Te=%g' % (k, Ti_Te))

subplot(1,2,1)
contour(XX, YY, real(iawDisp), colors='k', levels=[0], linewidth=1)
contour(XX, YY, imag(iawDisp), colors='m', levels=[0], linewidth=1)
fAbs = numpy.ma.masked_less(abs(iawDisp), 10.0)
pcolormesh(XX, YY, log(fAbs))
grid()
axis('tight')

subplot(1,2,2)
contour(XX, YY, real(iawKinDisp), colors='k', levels=[0], linewidth=1)
contour(XX, YY, imag(iawKinDisp), colors='m', levels=[0], linewidth=1)
fAbs = numpy.ma.masked_less(abs(iawKinDisp), 10.0)
pcolormesh(XX, YY, log(fAbs))
grid()
axis('tight')

show()
