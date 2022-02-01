from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(pre, i):
    print("Working on %d ..." % i)
    
    data = pg.GData("%s_ions_M0_%d.bp" % (pre, i))
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate()

    data = pg.GData("%s_ions_M1i_%d.bp" % (pre, i))
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, nux = dg.interpolate(0)
    XX, nuy = dg.interpolate(1)    
    ux = nux/q0
    uy = nuy/q0
    
    data = pg.GData("%s_ions_M2_%d.bp" % (pre, i))
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q2 = dg.interpolate()

    dx = 1/q2.shape[0]

    vt = (q2-q0*(ux**2+uy**2))*0.5
    return dx*vt.sum()

figure(1)
Ts = []
for i in range(101):
    Ts.append(plotFig("s4-oscc-E", i))

Tc = []    
for i in range(101):
    Tc.append(plotFig("../c4/c4-oscc-E", i))

tm = linspace(0, 100, 101)    
plot(tm, Ts, 'r-', label='SC')
plot(tm, Tc, 'k-', label='TD')
xlabel('Time')
ylabel(r'$n v_{th}^2$')
legend(loc='best')
grid()
savefig('c4-temp-hist.png', dpi=150)

show()

