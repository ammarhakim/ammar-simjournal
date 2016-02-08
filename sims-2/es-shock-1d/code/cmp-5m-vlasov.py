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

from optparse import OptionParser

# set command line options
parser = OptionParser()
parser.add_option('-k', '--kinetic', action = 'store',
                  dest = 'kinFileName',
                  help = 'Hdf5 file prefix from kinetic sim to plot')
parser.add_option('-f', '--fluid', action = 'store',
                  dest = 'fluFileName',
                  help = 'Hdf5 file prefix from fluid sim to plot')
parser.add_option('-o', '--outprefix', action = 'store',
                  dest = 'outPrefix',
                  help = 'Output prefix for figures')
parser.add_option('-m', '--mass-ratio', action = 'store',
                  dest = 'massRatio',
                  help = 'Ion/electron mass ratio')
parser.add_option('-n', '--frame', action = 'store',
                  dest = 'frame',
                  help = 'Frame to plot')

(options, args) = parser.parse_args()

kinFileName = options.kinFileName
fluFileName = options.fluFileName
outPrefix = ''
if options.outPrefix:
    outPrefix = options.outPrefix
frame = int(options.frame)
massRatio = float(options.massRatio)

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
d = gkedata.GkeData("%s_numDensityElc_%d.h5" % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, numElc = dg.project(0)

tm = d.time # time from

d = gkedata.GkeData("%s_q_%d.h5" % (fluFileName, frame) )
fv = gkedgbasis.GkeDgPolyOrder0Basis(d)
XcFV, numElcFV = fv.project(0)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(1)
plot(Xc, numElc, 'r-')
plot(XcFV, numElcFV, 'k-')
title('Electron number density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-elc-numDens_%d.png' % (outPrefix, frame) )

# ion number density
d = gkedata.GkeData("%s_numDensityIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, numIon = dg.project(0)

XcFV, numIonFV = fv.project(5)
XcFV = XcFV/XcFV[-1]*Xc[-1]
numIonFV = numIonFV/massRatio

figure(2)
plot(Xc, numIon, 'r-')
plot(XcFV, numIonFV, 'k-')
title('Ion number density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-ion-numDens_%d.png' % (outPrefix, frame) )

######

# electron momentum density
d = gkedata.GkeData("%s_momentumElc_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momElc = dg.project(0)

XcFV, momElcFV = fv.project(1)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(3)
plot(Xc, momElc, 'r-')
plot(XcFV, momElcFV/0.01, 'k-')
title('Electron momentum density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-elc-momDens_%d.png' % (outPrefix, frame) )

# ion momentum density
d = gkedata.GkeData("%s_momentumIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momIon = dg.project(0)

XcFV, momIonFV = fv.project(6)
XcFV = XcFV/XcFV[-1]*Xc[-1]
momIonFV = momIonFV/massRatio

figure(4)
plot(Xc, momIon, 'r-')
plot(XcFV, momIonFV/0.01, 'k-')
title('Ion momentum density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-ion-momDens_%d.png' % (outPrefix, frame) )

######

# electron energy density
d = gkedata.GkeData("%s_ptclEnergyElc_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erElc = dg.project(0)
erElc = 0.5*erElc

XcFV, erElcFV = fv.project(4)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(5)
plot(Xc, erElc, 'r-')
plot(XcFV, erElcFV/erElcFV[0]*erElc[0], 'k-')
title('Electron energy density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-elc-enrDens_%d.png' % (outPrefix, frame) )

# ion momentum density
d = gkedata.GkeData("%s_ptclEnergyIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erIon = dg.project(0)
erIon = 0.5*massRatio*erIon

XcFV, erIonFV = fv.project(9)
XcFV = XcFV/XcFV[-1]*Xc[-1]
erIonFV = erIonFV

figure(6)
plot(Xc, erIon, 'r-')
plot(XcFV, erIonFV/erIonFV[0]*erIon[0], 'k-')
title('Ion energy density at t=%g' % tm)
savefig('%s-fluid-kinetic-cmp-ion-enrDens_%d.png' % (outPrefix, frame) )

#show()
