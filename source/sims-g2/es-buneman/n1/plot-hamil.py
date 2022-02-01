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

    return phi-0*min(phi)

def getPhiAndDistf(fr):
    data = pg.GData("n1-es-buneman_elc_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XXv, q0 = dg.interpolate()

    return XXv, q0, getPhi(fr)

# field
XXv, fe, phi = getPhiAndDistf(86)
X, Y = meshgrid(XXv[0], XXv[1])

phi2, V2 = meshgrid(phi, XXv[1])
H = 0.5*V2**2 - phi2
    
fig, ax = subplots(1,1)
im = ax.pcolormesh(X, Y, transpose(fe[:,:,0]), vmin=0.0)
#ax.contour(X, Y, H, levels=[0], colors='w', linewidths=1)
ax.contour(X, Y, H, 40, colors='w', linewidths=1)

show()
