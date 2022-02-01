from pylab import *
import postgkyl as pg
import numpy

def _colorbar(obj, fig, ax, label=""):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(obj, cax=cax, label=label)

style.use('postgkyl.mplstyle')

for i in range(69):
    print("Working on %d ... " % i)

    fig = figure(1)

    # ions
    data = pg.GData("n1-Lz4-no-collisions_ion_GkM0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 1, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    nz2 = int(q0.shape[2]/2)
    q0m = numpy.ma.masked_where(q0[:,:,nz2,0]<0, q0[:,:,nz2,0])
    vmax = q0m.max()
    vmin = q0m.min()
    subplot(1,2,1)
    im = pcolormesh(X, Y, transpose(q0m))
    title("Ion")
    axis('image')

    # elc
    data = pg.GData("n1-Lz4-no-collisions_electron_GkM0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 1, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    nz2 = int(q0.shape[2]/2)
    q0m = numpy.ma.masked_where(q0[:,:,nz2,0]<0, q0[:,:,nz2,0])
    subplot(1,2,2)
    im = pcolormesh(X, Y, transpose(q0m), vmin=vmin, vmax=vmax)
    title("Electron")
    axis('image')    

    suptitle(r'Time %g $\mu s$' % (i*1.0))

    savefig("n1-Lz4-no-collisions_GkM0_%05d.png" % i, dpi=200)
    close()
    

