import gkedata
import gkedgbasis
from pylab import *

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
rcParams['xtick.major.size']           = 8 # default is 4
rcParams['xtick.major.width']          = 3 # default is 0.5
rcParams['ytick.major.size']           = 8 # default is 4
rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
rcParams['savefig.bbox']               = 'tight'

# Math/LaTex fonts:
# http://matplotlib.org/users/mathtext.html
# http://matplotlib.org/users/usetex.html
# Example: xlabel(r'$t \cdot l / V_{A,bc}$')
rcParams['mathtext.default'] = 'regular' # match the font used for regular text

# electron number density
d = gkedata.GkeData("../s1/s1-es-shock_numDensityElc_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, numElc = dg.project(0)

d = gkedata.GkeData("s2-5m-es-shock_q_100.h5")
fv = gkedgbasis.GkeDgPolyOrder0Basis(d)
XcFV, numElcFV = fv.project(0)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(1)
plot(Xc, numElc, 'r-')
plot(XcFV, numElcFV, 'k-')
title('Electron number density')
savefig('fluid-kinetic-cmp-elc-numDens.png')

# ion number density
d = gkedata.GkeData("../s1/s1-es-shock_numDensityIon_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, numIon = dg.project(0)

XcFV, numIonFV = fv.project(5)
XcFV = XcFV/XcFV[-1]*Xc[-1]
numIonFV = numIonFV/1836.2

figure(2)
plot(Xc, numIon, 'r-')
plot(XcFV, numIonFV, 'k-')
title('Ion number density')
savefig('fluid-kinetic-cmp-ion-numDens.png')

######

# electron momentum density
d = gkedata.GkeData("../s1/s1-es-shock_momentumElc_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momElc = dg.project(0)

d = gkedata.GkeData("s2-5m-es-shock_q_100.h5")
fv = gkedgbasis.GkeDgPolyOrder0Basis(d)
XcFV, momElcFV = fv.project(1)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(3)
plot(Xc, momElc, 'r-')
plot(XcFV, momElcFV/0.01, 'k-')
title('Electron momentum density')
savefig('fluid-kinetic-cmp-elc-momDens.png')

# ion momentum density
d = gkedata.GkeData("../s1/s1-es-shock_momentumIon_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momIon = dg.project(0)

XcFV, momIonFV = fv.project(6)
XcFV = XcFV/XcFV[-1]*Xc[-1]
momIonFV = momIonFV/1836.2

figure(4)
plot(Xc, momIon, 'r-')
plot(XcFV, momIonFV/0.01, 'k-')
title('Ion momentum density')
savefig('fluid-kinetic-cmp-ion-momDens.png')

######

# electron energy density
d = gkedata.GkeData("../s1/s1-es-shock_ptclEnergyElc_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erElc = dg.project(0)
erElc = 0.5*erElc

d = gkedata.GkeData("s2-5m-es-shock_q_100.h5")
fv = gkedgbasis.GkeDgPolyOrder0Basis(d)
XcFV, erElcFV = fv.project(4)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(5)
plot(Xc, erElc, 'r-')
plot(XcFV, erElcFV/erElcFV[0]*erElc[0], 'k-')
title('Electron energy density')
savefig('fluid-kinetic-cmp-elc-enrDens.png')

# ion momentum density
d = gkedata.GkeData("../s1/s1-es-shock_ptclEnergyIon_100.h5")
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erIon = dg.project(0)
erIon = 0.5*1836.2*erIon

XcFV, erIonFV = fv.project(9)
XcFV = XcFV/XcFV[-1]*Xc[-1]
erIonFV = erIonFV

figure(6)
plot(Xc, erIon, 'r-')
plot(XcFV, erIonFV/erIonFV[0]*erIon[0], 'k-')
title('Ion energy density')
savefig('fluid-kinetic-cmp-ion-enrDens.png')

show()
