from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(i,fr):
    print("Working on %d ..." % i)
    figure(i)
    data = pg.GData("c6-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    nx2 = int(q.shape[0]/2)
    nvx2 = int(q.shape[1]/2)    

    subplot(2,1,1)
    pcolormesh(XX[1], XX[2], transpose(q[nx2,:,:,0]))
    grid()
    axis('image')

    subplot(2,1,2)
    pcolormesh(XX[0], XX[2], transpose(q[:,nvx2,:,0]))
    grid()

    savefig("oscc-E-mov-%05d.png" % i)
    close()    

for i in range(101):
    plotFig(i,i)
    
show()

