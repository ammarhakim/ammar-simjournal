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

X, s = calcMeanSpect('s5-euler-kh_chi-spect', 190, 200)
loglog(X, s/s[0])
grid()

xx = linspace(13, 100, 20)
ee = 0.75*(xx/xx[0])**(-1.0)
loglog(xx, ee, '--k')
text(50,0.5,r'$k^{-1}$')

xlabel(r'$k/2\pi$')
ylabel(r'$\chi$')
savefig('chi-spect.png')

show()
    
