from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i)
    data = pg.GData("c4-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    qSum = sum(q, axis=0)
    
    pcolormesh(XX[1], XX[2], transpose(qSum[:,:,0]))
    grid()
    axis('image')
    savefig("oscc-E-%05d.png" % i, dpi=150)
    close()

for i in range(101):
    plotFig(i,i)
    
show()

