from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

E0 = 0.5

def wR(t):
    wx = E0/2*(t*cos(t)+sin(t))
    wy = -E0/2*t*sin(t)
    return wx, wy

T = linspace(0, 20, 500)
wx, wy = wR(T)

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i)
    data = pg.GData("s3-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()

    pcolormesh(XX[1], XX[2], transpose(q[3,:,:,0]))
    #colorbar()
    plot(wx, wy, 'w', linewidth=1)
    grid()
    axis('image')
    savefig("oscc-E-%05d.png" % i, dpi=150)
    close()

for i in range(21):
    plotFig(i,i)
    
show()

