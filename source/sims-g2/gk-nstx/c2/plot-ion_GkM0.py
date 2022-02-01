from pylab import *
import postgkyl as pg
import numpy

def _colorbar(obj, fig, ax, label=""):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(obj, cax=cax, label=label)

style.use('postgkyl.mplstyle')

for i in range(100):
    print("Working on %d ... " % i)

    fig, ax = subplots(1,2)

    # ions
    data = pg.GData("c2-Lz4-lbo-collisions_ion_GkM0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 1, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    nz2 = int(q0.shape[2]/2)
    q0m = numpy.ma.masked_where(q0[:,:,nz2,0]<0, q0[:,:,nz2,0])
    vmax = q0m.max()
    vmin = q0m.min()
    im = ax[0].pcolormesh(X, Y, transpose(q0m))
    _colorbar(im, fig, ax[0])
    ax[0].axis('image')

    # elc
    data = pg.GData("c2-Lz4-lbo-collisions_electron_GkM0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 1, "ms")
    XX, q0 = dg.interpolate()
    X, Y = meshgrid(XX[0], XX[1])
    nz2 = int(q0.shape[2]/2)
    q0m = numpy.ma.masked_where(q0[:,:,nz2,0]<0, q0[:,:,nz2,0])
    im = ax[1].pcolormesh(X, Y, transpose(q0m), vmin=vmin, vmax=vmax)
    ax[1].set_yticks([])
    
    _colorbar(im, fig, ax[1])
    ax[1].axis('image')

    suptitle(r'Time %g $\mu s$' % (i*1.0))

    savefig("c2-Lz4-lbo-collisions_GkM0_%05d.png" % i, dpi=200)
    close()
    

