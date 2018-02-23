from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(i, c, lbl):
    print("Working on %d ..." % i)
    
    data = pg.GData("c5-oscc-E_ions_M0_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate()

    data = pg.GData("c5-oscc-E_ions_M1i_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, nux = dg.interpolate(0)
    XX, nuy = dg.interpolate(1)    
    ux = nux/q0
    uy = nuy/q0
    
    data = pg.GData("c5-oscc-E_ions_M2_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q2 = dg.interpolate()

    dx = 1/q2.shape[0]

    vt = (q2-q0*(ux**2+uy**2))*0.5
    plot(XX[0], vt, c, label=lbl)
    plot([XX[0][0], XX[0][-1]], [dx*vt.sum(), dx*vt.sum()], "%s--" % c)

figure(1)
plotFig(0, "k", "$t=0$")
plotFig(25, "r", "$t=25$")
plotFig(50, "m", "$t=50$")
plotFig(100, "b", "$t=100$")
legend(loc='best')
grid()
xlabel("$X$")
ylabel("$n v_{th}^2$")

savefig('c5-temp-cmp-1d.png', dpi=150)

show()

