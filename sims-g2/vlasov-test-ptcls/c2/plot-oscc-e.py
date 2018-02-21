from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

E0 = 1.0
w = 0.5

def wNR(t):
    wx = E0/(1-w**2)*(sin(t)-w*sin(w*t))
    wy = E0/(1-w**2)*(cos(t)-cos(w*t))
    return wx, wy

T = linspace(0, 100, 5000)
wx, wy = wNR(T)

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i)
    data = pg.GData("c2-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()

    pcolormesh(XX[1], XX[2], transpose(q[3,:,:,0]))
    #colorbar()
    plot(wx, wy, 'w', linewidth=1)
    axis('image')
    savefig("oscc-E-%05d.png" % i, dpi=150)
    close()

for i in range(101):
    plotFig(i,i)
    
show()

