from pylab import *
import postgkyl

style.use('postgkyl.mplstyle')

def calcMeanSpect(fn, lo, up):
    d = postgkyl.GData("%s_00%d.bp" % (fn, lo))
    s = d.getValues()
    g = d.getGrid()
    dx = g[0][1]-g[0][0]
    gc = linspace(g[0][0]+0.5*dx, g[0][-1]-0.5*dx, g[0].shape[0]-1)

    for i in range(lo+1, up+1):
        d = postgkyl.GData("%s_00%d.bp" % (fn, i))
        s = s + d.getValues()

    return gc, sqrt(s/(up-lo+1))

ng1, n1 = calcMeanSpect("../n1/n1-euler-kh_ke-spect", 190, 200)
sg1, s1 = calcMeanSpect("../s1/s1-euler-kh_ke-spect", 190, 200)

loglog(ng1, n1/n1[0], label=r'Random')
loglog(sg1, s1/s1[0], label=r'1-Mode')
legend(loc='best')
grid()

title('Kinetic energy spectrum')
xlabel(r'$k/2\pi$')
ylabel(r'$KE_k/KE_{k_{min}}$')
savefig('ke-spect-n1-s1-cmp.png', dpi=150)

show()
