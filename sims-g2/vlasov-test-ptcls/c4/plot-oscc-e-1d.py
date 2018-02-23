from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(i,fr):
    print("Working on %d ..." % i)
    data = pg.GData("c4-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    nx2 = int(q.shape[0]/2)
    nvx2 = int(q.shape[1]/2)    

    return q[nx2,nvx2,:,0]


nFrame = 21

q1 = zeros((nFrame+1,24*3), float)
for i in range(nFrame+1):
    q1[i,:] = plotFig(i,i)

contour(q1, colors='k')
#pcolormesh(q1)
show()

