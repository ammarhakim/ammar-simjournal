from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

w = 0.4567
tper = 2*pi/w

def getData(i):
    #print("Working on %d ..." % i)
    data = pg.GData("c4-oscc-E_ions_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    nx, nvx, nvy = q.shape[0], q.shape[1], q.shape[2]
    
    qMid = q[int(nx/2),:,int(nvy/2)]

    tm = data.time
    tm = tm - tper*floor(tm/tper)

    return tm, qMid

nFrame = 100
qFull = zeros((nFrame+1,72), float)
T = zeros((nFrame+1,), float)

for i in range(nFrame+1):
    tm, qt = getData(i)
    T[i] = tm
    
    qFull[i,:] = squeeze(qt)

# sort to enforce periodicity in time
idx = T.argsort()
Ts = T[idx]

qS = 0.0*qFull
for i in range(idx.shape[0]):
    qS[i,:] = qFull[idx[i],:]

def maxwellian2D(n, vx, vy, ux, uy, vth):
   v2 = (vx - ux)**2 + (vy - uy)**2
   return n/(2*pi*vth**2)*exp(-v2/(2*vth**2))

Vx = linspace(-6, 6, 72)
contour(Ts, Vx, transpose(qS), colors='k', linestyles='dotted')
savefig("c4-pc.png", dpi=150)
show()
