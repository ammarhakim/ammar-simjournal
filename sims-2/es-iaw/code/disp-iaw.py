from scipy import special
import math
from pylab import *
import numpy.ma

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
Ti_Te = 0.5
gasGamma = 1.0
k = 0.5

Ti = mi*vthi**2
Te = Ti/Ti_Te

me = mi*mu

wpe = wpi/math.sqrt(mu)
cse = sqrt(gasGamma)*vthi/math.sqrt(Ti_Te)/math.sqrt(mu)
csi = sqrt(3)*vthi

ionSound = 1.040707021951675
print("Plasma frequency = %g. Electron sound speed = %g" % (wpe, cse))
print("Ion sound speed %g" % ionSound)

def iawEps(w, k):
    return 1-wpi**2/(2*vthi**2*k**2)*derivPlasmaDisp(w/(math.sqrt(2)*vthi*k)) - wpe**2/(w**2-cse**2*k**2)

def iawFluid(w, k):
    return 1-wpi**2/(w**2-csi**2*k**2) - wpe**2/(w**2-cse**2*k**2)

X = linspace(-5.0, 5.0, 400)
Y = linspace(-2.0, 0.1, 400)
XX, YY = meshgrid(X, Y)
ZZ = XX+1j*YY

figure(1)
iawDisp = iawEps(ZZ, k)
iawFluidDisp = iawFluid(ZZ, k)

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

contour(XX, YY, real(iawDisp), colors='k', levels=[0])
contour(XX, YY, imag(iawDisp), colors='m', levels=[0])
fAbs = numpy.ma.masked_less(abs(iawDisp), 10.0)
pcolormesh(XX, YY, log(fAbs))
plot([-ionSound, -ionSound], gca().get_ylim(), 'r--')
plot([ionSound, ionSound], gca().get_ylim(), 'r--')

grid()
axis('tight')

show()
