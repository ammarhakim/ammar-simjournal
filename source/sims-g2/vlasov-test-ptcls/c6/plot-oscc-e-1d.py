from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(i,fr):
    print("Working on %d ..." % i)
    data = pg.GData("c6-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    qsum = sum(q,axis=0)

    nvx2 = int(q.shape[1]/2)

    return qsum[nvx2,:].squeeze()


nFrame = 100

q1 = zeros((nFrame+1,24*3), float)
for i in range(nFrame+1):
    q1[i,:] = plotFig(i,i)
    plot(q1[i,:])

figure(2)    
contour(q1, colors='k')
#pcolormesh(q1)
#colorbar()
show()

