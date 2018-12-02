from pylab import *
import postgkyl as pg
import numpy

style.use('postgkyl.mplstyle')

for i in range(0, 151):
    print("Working on %d ... " % i)

    fig, ax = subplots(2,1)

    # ions
    data = pg.GData("c1-buneman_ion_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    im = ax[0].pcolormesh(X, Y, transpose(q0[:,:,0]), vmin=0.0)

    # elc
    data = pg.GData("c1-buneman_elc_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    im = ax[1].pcolormesh(X, Y, transpose(q0[:,:,0]), vmin=0.0)
    
    suptitle(r'Time %g' % (i*1.0))

    savefig("c1-buneman_%05d.png" % i, dpi=200)
    close()

