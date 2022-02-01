import numpy
from pylab import *
import tables

rc('text', usetex=True)

def findXloc(X, fx, val):
    cross = []
    for i in range(X.shape[0]-1):
        if val>fx[i] and val<fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))
        elif val<fx[i] and val>fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))            
    return cross

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*elcMass*n0)
ionCycl = ionCharge*B0/ionMass

timeSlice = []
extent = []

for i in range(15,41):
    print ('Working on %d ...' % i)
    fh = tables.openFile ("s279-gem-tenmom_q_%d.h5" % i)
    tm = fh.root.timeData._v_attrs['vsTime']*ionCycl
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]

    X = linspace(-LX/2, LX/2, nx)
    Y = linspace(-LY/2, LY/2, ny)

    ux = q[:,:,1]/q[:,:,0]
    ux1 = ux[:,ny/2]
    minix = ux1.argmin()
    maxix = ux1.argmax()
    minx = X[minix]
    maxx = X[maxix]

    # store time and current sheet extent
    timeSlice.append(fh.root.timeData._v_attrs['vsTime']*ionCycl)
    extent.append(maxx-minx)

timeSliceA = numpy.array(timeSlice)
extentA = numpy.array(extent)

plot(timeSliceA, extentA, '-ro')
title('10M Current sheet extension v/s time')
xlabel('Time ($t\Omega_{ci}$)')
ylabel('Extension ($d_i$)')
savefig('s279-10m-current-extension-vs-t.png')

show()
