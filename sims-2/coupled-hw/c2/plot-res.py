from pylab import *
import postgkyl as pg

def colorbar_adj(obj, mode=1, redraw=False, _fig_=None, _ax_=None, aspect=None):
    '''
    Add a colorbar adjacent to obj, with a matching height
    For use of aspect, see http://matplotlib.org/api/axes_api.html#matplotlib.axes.Axes.set_aspect ; E.g., to fill the rectangle, try "auto"
    '''
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    if mode == 1:
        _fig_ = obj.figure; _ax_ = obj.axes
    elif mode == 2: # assume obj is in the current figure/axis instance
        _fig_ = plt.gcf(); _ax_ = plt.gca()
    _divider_ = make_axes_locatable(_ax_)
    _cax_ = _divider_.append_axes("bottom", size="5%", pad=0.3)
    _cbar_ =  _fig_.colorbar(obj, cax=_cax_, orientation='horizontal')
    if aspect != None:
        _ax_.set_aspect(aspect)
    if redraw:
        _fig_.canvas.draw()
    return _cbar_

# Early time
d = pg.GData("c2-hw_chiFull_5.h5")
dgC = pg.data.interp.GInterpNodalSerendipity(d, 2)
X, chiC = dgC.project(0)

d = pg.GData("../f2/f2-hw_chi_5.h5")
dgS = pg.data.interp.GInterpNodalSerendipity(d, 2)
X, chiS = dgS.project(0)

figure(1)
suptitle('Vorticity: Global (left), Coupled (right). T=5')
fig = subplot(1,2,1)
im = pcolormesh(X[0], X[1], transpose(chiS), cmap='inferno')
axis('image')
colorbar_adj(im)

subplot(1,2,2)
im = pcolormesh(X[0], X[1], transpose(chiC), cmap='inferno')
plot([0, 0], gca().get_ylim(), 'w--')
plot([5, 5], gca().get_ylim(), 'w--')
axis('image')
colorbar_adj(im)

savefig('c2-coupled-chi-cmp-t5.png', dpi=300)

# Late time
d = pg.GData("c2-hw_chiFull_10.h5")
dgC = pg.data.interp.GInterpNodalSerendipity(d, 2)
X, chiC = dgC.project(0)

d = pg.GData("../f2/f2-hw_chi_10.h5")
dgS = pg.data.interp.GInterpNodalSerendipity(d, 2)
X, chiS = dgS.project(0)

figure(2)
suptitle('Vorticity: Global (left), Coupled (right). T=10')
fig = subplot(1,2,1)
im = pcolormesh(X[0], X[1], transpose(chiS), cmap='inferno')
axis('image')
colorbar_adj(im)

subplot(1,2,2)
im = pcolormesh(X[0], X[1], transpose(chiC), cmap='inferno')
plot([0, 0], gca().get_ylim(), 'w--')
plot([5, 5], gca().get_ylim(), 'w--')
axis('image')
colorbar_adj(im)

savefig('c2-coupled-chi-cmp-t10.png', dpi=300)

show()


