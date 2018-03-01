from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

def plotFig(nr,nc,i,fr):
    print("Working on %d ..." % i)
    
    data = pg.GData("s5-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    qSum = sum(q, axis=0)

    subplot(nr,nc,i)
    pcolormesh(XX[1], XX[2], transpose(qSum[:,:,0]))
    grid()
    axis('image')

figure(1)
plotFig(2,2,1,0)
ylabel("$v_y$")

plotFig(2,2,2,25)

plotFig(2,2,3,50)
xlabel("$v_x$")
ylabel("$v_y$")

plotFig(2,2,4,100)
xlabel("$v_x$")

savefig('s5-vxvy-cmp.png', dpi=150)

def plotFig1D(i, lbl):
    print("Working on %d ..." % i)
    
    data = pg.GData("s5-oscc-E_ions_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()
    numCells = q.shape
    qSum = sum(q, axis=0)*(2.0*pi)/numCells[0]

    plot(XX[1], qSum[:,int(numCells[2]/2),0], label=lbl)
    grid()

figure(2)
plotFig1D(0, "$t=0$")
plotFig1D(25, "$t=25$")
plotFig1D(50, "$t=50$")
plotFig1D(100, "$t=100$")
legend(loc="best")
xlabel('$v_x$')
ylabel('$f(v_x,v_y=0,t)$')

savefig('s5-vxvy-cmp-1d.png', dpi=150)

show()

