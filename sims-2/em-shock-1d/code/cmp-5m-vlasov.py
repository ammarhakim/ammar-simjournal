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
rcParams['lines.linewidth']            = 1 #2
rcParams['font.size']                  = 12
rcParams['xtick.major.size']           = 4 #8 # default is 4
rcParams['xtick.major.width']          = 0.5 # default is 0.5
rcParams['ytick.major.size']           = 4 # default is 4
rcParams['ytick.major.width']          = 0.5 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.25 #0.125
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

figure(1)

d = gkedata.GkeData("%s_q_%d.h5" % (fluFileName, frame) )
fv = gkedgbasis.GkeDgPolyOrder0Basis(d)
XcFV, numElcFV = fv.project(0)
XcFV = XcFV/XcFV[-1]*Xc[-1]

subplot(3,2,1)
title('Electrons')
plot(Xc, numElc, 'r-')
plot(XcFV, numElcFV, 'k-')
ylabel('Density')
#title('Electron number density at t=%g' % tm)

# ion number density
d = gkedata.GkeData("%s_numDensityIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, numIon = dg.project(0)

XcFV, numIonFV = fv.project(5)
XcFV = XcFV/XcFV[-1]*Xc[-1]
numIonFV = numIonFV/massRatio

ax = subplot(3,2,2)
title('Ions')
plot(Xc, numIon, 'r-')
plot(XcFV, numIonFV, 'k-')
#title('Ion number density at t=%g' % tm)

######

# electron momentum density
d = gkedata.GkeData("%s_momentumElc_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momElc = dg.project(0)

XcFV, momElcFV = fv.project(1)
XcFV = XcFV/XcFV[-1]*Xc[-1]

subplot(3,2,3)
plot(Xc, momElc, 'r-')
plot(XcFV, momElcFV, 'k-')
ylabel('Momentum')
#title('Electron momentum density at t=%g' % tm)

# ion momentum density
d = gkedata.GkeData("%s_momentumIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, momIon = dg.project(0)

XcFV, momIonFV = fv.project(6)
XcFV = XcFV/XcFV[-1]*Xc[-1]
momIonFV = momIonFV/massRatio

subplot(3,2,4)
plot(Xc, momIon, 'r-')
plot(XcFV, momIonFV, 'k-')
#title('Ion momentum density at t=%g' % tm)

######

# electron energy density
d = gkedata.GkeData("%s_ptclEnergyElc_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erElc = dg.project(0)
erElc = 0.5*erElc

XcFV, erElcFV = fv.project(4)
XcFV = XcFV/XcFV[-1]*Xc[-1]

subplot(3,2,5)
plot(Xc, erElc, 'r-')
plot(XcFV, erElcFV/erElcFV[0]*erElc[0], 'k-')
ylabel('Energy')
#title('Electron energy density at t=%g' % tm)

# ion momentum density
d = gkedata.GkeData("%s_ptclEnergyIon_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, erIon = dg.project(0)
erIon = 0.5*massRatio*erIon

XcFV, erIonFV = fv.project(9)
XcFV = XcFV/XcFV[-1]*Xc[-1]
erIonFV = erIonFV

subplot(3,2,6)
plot(Xc, erIon, 'r-')
plot(XcFV, erIonFV/erIonFV[0]*erIon[0], 'k-')
#title('Ion energy density at t=%g' % tm)

suptitle("Kinetic (red) and Fluid (black) solution at t=%g" % tm)

savefig('%s-cmp-kin-flu_%d.png' % (outPrefix, frame), dpi=300)

d = gkedata.GkeData("%s_em_%d.h5"  % (kinFileName, frame) )
dg = gkedgbasis.GkeDgLobatto1DPolyOrder2Basis(d)
Xc, Ex = dg.project(0)
Xc, Bz = dg.project(5)

XcFV, exFV = fv.project(10)
XcFV = XcFV/XcFV[-1]*Xc[-1]

XcFV, bzFV = fv.project(15)
XcFV = XcFV/XcFV[-1]*Xc[-1]

figure(2)
plot(Xc, Ex, 'r-')
plot(XcFV, exFV, 'k-')
savefig('%s-cmp-kin-flu-ex_%d.png' % (outPrefix, frame), dpi=300)

figure(3)
plot(Xc, Bz, 'r-')
plot(XcFV, bzFV, 'k-')
savefig('%s-cmp-kin-flu-bz_%d.png' % (outPrefix, frame), dpi=300)

#show()
