from pylab import *
import postgkyl as pg
import numpy
import numpy as np

style.use('postgkyl.mplstyle')

for i in range(0, 39):
    print("Working on %d ... " % i)

    fig = figure(1)

    data = pg.GData("r2-vlasov-rt-mass100_elc_M0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate(0)
    X, Y = meshgrid(XX[0], XX[1])
    im = pcolormesh(X, Y, transpose(q0[:,:,0]))

    title(r"$\Omega t = {:0.2f}$".format(60*data.time/60000))
    axis('image')

    savefig("r2-vlasov-rt-mass100_elc_M0_%05d.png" % i, dpi=200)
    close()
    

