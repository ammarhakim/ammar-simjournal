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

g3, s3 = calcMeanSpect("../s3/s3-euler-kh_ke-spect", 190, 200)
g5, s5 = calcMeanSpect("../s5/s5-euler-kh_ke-spect", 190, 200)

loglog(g3, s3/s3[0], label=r'$2000\times 2000$')
loglog(g5, s5/s5[0], label=r'$4096\times 4096$')
legend(loc='best')
grid()

xx = linspace(5, 500, 20)
ee = 1/10.0*(xx/xx[0])**(-2.0)
loglog(xx, ee, '--r')
text(50,5e-3,r'$k^{-2}$')

xx2 = linspace(500, 2000, 20)
ee2 = 6.5e-6*(xx2/xx2[0])**(-3.5)
loglog(xx2, ee2, '--m')
text(1000,1.5e-6,r'$k^{-3.5}$')

xx3 = linspace(200, 1000, 20)
ee3 = 6.5e-6*(xx3/xx3[0])**(-3.5)
loglog(xx3, ee3, '--m')

title('Kinetic energy spectrum')
xlabel(r'$k/2\pi$')
ylabel(r'$KE_k/KE_{k_{min}}$')
savefig('ke-spect-s3-s5-cmp.png', dpi=150)

show()
