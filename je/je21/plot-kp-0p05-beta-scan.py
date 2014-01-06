from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import numpy
import math

# customization for figure
rcParams['lines.linewidth']            = 1.5
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

fig, ax1 = subplots()

# plot frequencies
dat = loadtxt('kp-0p05-beta-scan')
freq = 2*pi/(2*dat[:,1])
ax1.plot(dat[:,0], freq/(math.sqrt(2)*1.0*0.5), 'mo', markersize=10)
#ax1.semilogx(dat[:,0], freq/(math.sqrt(2)*1.0*0.5), 'mo', markersize=10)
ax1.set_xlabel(r'$\hat{\beta}$')
datEx = loadtxt('KAW-kp-0p05-rates.txt')
ax1.plot(datEx[:,0], datEx[:,1], 'm-', label='Exact')
#ax1.semilogx(datEx[:,0], datEx[:,1], 'm-', label='Exact')
ax1.set_ylabel('Normalized Frequency ($\Omega/\sqrt{2} k_{\parallel}$)', color='m')

ax2 = ax1.twinx()
ax2.plot(dat[:,0], 0.5*dat[:,2]/(math.sqrt(2)*1.0*0.5), 'go', markersize=10)
#ax2.semilogx(dat[:,0], 0.5*dat[:,2]/(math.sqrt(2)*1.0*0.5), 'go', markersize=10)
ax2.plot(datEx[:,0], -datEx[:,2], 'g-', label='Exact')
#ax2.semilogx(datEx[:,0], -datEx[:,2], 'g-', label='Exact')
ax2.set_ylabel('Normalized Damping ($\gamma/\sqrt{2} k_{\parallel}$)', color='g')
title('$k_\perp^2=0.05$')

pylab.axis('tight')
savefig('freq-damp-shear-alf-kp-0p05-beta-scan.png', bbox_inches='tight')
show()
close()
