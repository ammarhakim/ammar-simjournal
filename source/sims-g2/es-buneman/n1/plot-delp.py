from pylab import *
import postgkyl as pg
import numpy

style.use('postgkyl.mplstyle')

def getPhi(fr):
    data = pg.GData("n1-es-buneman_field_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, E = dg.interpolate()
    X1 = XX[0]
    dx = X1[1]-X1[0]
    
    phi = numpy.zeros((E.shape[0]+1,1), numpy.float)
    for i in range(1,phi.shape[0]):
        phi[i] = phi[i-1]-dx*E[i-1]

    return phi

v = []
for i in range(0, 151):
    p = getPhi(i)
    v.append(100*(p[-1,0]-p[0,0])/abs(p).max())

plot(v)
xlabel('Time')
ylabel('% error in phi[Left]-phi[Right]')
grid()
savefig('error-in-phi.png')
show()
