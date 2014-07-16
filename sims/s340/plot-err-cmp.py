import tables
import pylab
import math
import numpy
from matplotlib import rcParams
import matplotlib.pyplot as plt

#pylab.rc('text', usetex=True)
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



from pylab import *

solL2 = 3.1399614 # this is L2 norm of initial condition
baseList = ["../s337/s337-modal-dg-diffuse", 
            "../s339/s339-modal-dg-diffuse", "../s340/s340-modal-dg-diffuse"]
titleStr = ["LDG-A", "LDG-S", "RDG"]
colStr = ['-g', '-k', '-r']
ordStr = ['2', '2', '4']

count = 0
figure(1)
for baseName in baseList:
    dat = loadtxt("%s_l2Norm.txt" % baseName)
    T = dat[:,0]
    val = dat[:,1]
    semilogy(T, sqrt(val/solL2), colStr[count], label=titleStr[count])
    count = count+1
xlabel('Time')
ylabel('RMS Error')
legend(loc='best')

savefig('s337-s338-s339-s340-rel-rms-err.pdf')

show()


