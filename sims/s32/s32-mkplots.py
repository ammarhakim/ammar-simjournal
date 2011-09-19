import numpy
import pylab
import tables

pylab.rc('text', usetex=True)
pylab.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,
                      wspace=0.2, hspace=0.2)

for m in range(2):
    # open data file from Garcia & Siewert paper
    gsDatNm = "gs-radiances-m%d.csv" % m
    gsDat = numpy.loadtxt(gsDatNm, delimiter=",")
    gsMu = gsDat[:,0]

    # now open HDF5 file
    lcFh = tables.openFile("s32-rte-slab_sol.h5")
    mu = lcFh.root.mu.read()

    muM = mu*0.0
    for i in range(mu.shape[0]):
        muM[i] = -mu[mu.shape[0]-i-1]

    muExtended = numpy.zeros( (2*mu.shape[0], ), numpy.float )
    muExtended[0:mu.shape[0]] = muM
    muExtended[mu.shape[0]:] = mu

    radiances = numpy.zeros( (2*mu.shape[0], ), numpy.float )

    factor = 1.0
    if (m==0):
        factor = 2.0

    #pylab.figure(m)

    titleList = [r'$\tau=0$', r'$\tau=\tau_0/20$', r'$\tau=\tau_0/10$', \
                 r'$\tau=\tau_0/5$', r'\tau=\tau_0/2$', r'$\tau=3\tau_0/4$', r'$\tau=\tau_0$']
    for d in range(7):
        up = lcFh.root._v_children["upward_radiance_%d" % d].read()
        upM = up[m,:]*0.0
        for i in range(upM.shape[0]):
            upM[i] = up[m,upM.shape[0]-i-1]
        down = lcFh.root._v_children["downward_radiance_%d" % d].read()

        radiances[0:mu.shape[0]] = upM
        radiances[mu.shape[0]:] = down[m,:]

        pylab.figure(7*m+d)
        #pylab.subplot(7, 1, d)
        pylab.plot(gsMu, factor*gsDat[:,d+1], 'ro')
        pylab.plot(muExtended, radiances, '-k')
        pylab.xlabel(r"Cosine of polar angle ($\mu$)")
        pylab.ylabel("Radiance")
        pylab.title(titleList[d])
        pylab.savefig("s32-rte-slab-m%d-tau%d.png" % (m, d))
        

    
