from pylab import *
import postgkyl as pg

style.use('postgkyl.mplstyle')

w = 0.4567
tper = 2*pi/w

def mkPoincare(T, v):
    t = T[1:]
    x = v[1:,0]
    vx = v[1:,1]

    vp = []
    tp = []
    for i in range(x.shape[0]-1):
        if (x[i]<0 and x[i+1]>0) or (x[i]>0 and x[i+1]<0):
            tm = t[i]
            tm = tm - tper*floor(tm/tper)
            tp.append(tm)
            vp.append(vx[i])
    return tp, vp

def plotP(p, fName):
    for i in range(1,11):
        print("Woroking on %s_%d.bp ... " % (fName, i))
        data = pg.GData("%s_%d.bp" % (fName, i))
        T = data.peakGrid()
        v = data.peakValues()
        tp, vp = mkPoincare(T[0], v)
        p.scatter(2*pi*array(tp)/tper, array(vp), 0.1, marker='.')

fig, ax = subplots(2,2)

plotP(ax[0,0], "ptcl-tracing-E0_E0_2_ptcl")
ax[0,0].set_title('$E_0=0.2$')
ax[0,0].set_xticks([])

plotP(ax[1,0], "ptcl-tracing-E0_E0_4_ptcl")
ax[1,0].set_title('$E_0=0.4$')

plotP(ax[0,1], "ptcl-tracing-E0_E0_6_ptcl")
ax[0,1].set_title('$E_0=0.6$')
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])

plotP(ax[1,1], "ptcl-tracing-E0_E0_95_ptcl")
ax[1,1].set_title('$E_0=0.95$')
ax[1,1].set_yticks([])
    
savefig('ptcl-poincare-E0.png', dpi=300)
show()    
