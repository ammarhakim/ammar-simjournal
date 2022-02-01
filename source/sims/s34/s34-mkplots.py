import numpy
import pylab
import tables
import math
import matplotlib.transforms as mtransforms

pylab.rc('text', usetex=True)

# open HDF5 file
lcFh = tables.openFile("s34-rte-slab_sol.h5")
mu = lcFh.root.mu.read()

muM = mu*0.0
for i in range(mu.shape[0]):
    muM[i] = -mu[mu.shape[0]-i-1]

muExtended = numpy.zeros( (2*mu.shape[0], ), numpy.float )
muExtended[0:mu.shape[0]] = muM
muExtended[mu.shape[0]:] = mu

def computeRadiance(rads, phi):
    nm, nmu = rads.shape[0], rads.shape[1]
    totalRad = numpy.zeros( (nmu, ), numpy.float )
    for m in range(nm):
        factor = 1.0
        if m == 0: factor = 0.5
        totalRad = totalRad + factor*rads[m,:]*math.cos(m*phi)
    return totalRad

titleList = [r'$\tau=0$', r'$\tau=\tau_0/20$', r'$\tau=\tau_0/10$', \
             r'$\tau=\tau_0/5$', r'\tau=\tau_0/2$', r'$\tau=3\tau_0/4$', r'$\tau=\tau_0$']

phiVals = [0.0, math.pi/2]
fileNms = ["gs-radiances-phi0.csv", "gs-radiances-phiPi2.csv"]

count = 0
for phi in phiVals:

    gsDat = numpy.loadtxt(fileNms[count], delimiter=",")
    gsMu = gsDat[:,0]
    
    fig = pylab.figure(count)
    fig.subplots_adjust(hspace=1.0)

    radiances = numpy.zeros( (2*mu.shape[0], ), numpy.float )
    for d in range(7):
        up = lcFh.root._v_children["upward_radiance_%d" % d].read()
        down = lcFh.root._v_children["downward_radiance_%d" % d].read()

        totalUp = computeRadiance(up, phi)
        totalDown = computeRadiance(down, phi)

        totalUpM = totalUp*0.0
        for i in range(totalUp.shape[0]):
            totalUpM[i] = totalUp[totalUp.shape[0]-i-1]

        radiances[0:totalUp.shape[0]] = totalUpM
        radiances[totalUp.shape[0]:] = totalDown

        ax = pylab.subplot(7, 1, d+1)
        pylab.plot(gsMu, gsDat[:,d+1], 'ro')
        pylab.plot(muExtended, radiances, '-k')
        if d < 6:
            ax.set_xticklabels([""]) # zap labels from X axis
        ylims = ax.get_ylim()
        ax.set_yticks([ylims[0], 0.5*(ylims[0]+ylims[1]), ylims[1]])
        pylab.title(titleList[d])

    pylab.savefig("s34-rte-slab-%s.png" % fileNms[count])
    count = count + 1
    

    
    
