from pylab import *
import math

style.use('../code/postgkyl.mplstyle')

def readError(fname):
    dat = loadtxt(fname)
    N = dat[:,0]*1.0
    err = dat[:,1]
    err1 = err-err[-1]
    return N[:-1], err1[:-1]

N1, e1 = readError("../a4/error-dx.txt")
N2, e2 = readError("../a8/error-dx.txt")
N3, e3 = readError("../a12/error-dx.txt")

Nn = linspace(8,32,100)
a = 0.01**(1/4.0)*Nn[0]
f = (a/Nn)**4
loglog(N1, e1, 'r-', label='p=1')
loglog(Nn, f, 'r--', linewidth=1.0)
text(16,1.3e-3, r"$\Delta x^4$", color='r')

Nn = linspace(4,32,100)
a = 1.e-2**(1/5.0)*Nn[0]
f = (a/Nn)**5
loglog(N2, e2, 'k-', label='p=2')
loglog(Nn, f, 'k--', linewidth=1.0)
text(10,1.0e-4, r"$\Delta x^5$", color='k')

Nn = linspace(3,16,100)
a = 2.e-4**(1/6.0)*Nn[0]
f = (a/Nn)**6
loglog(N3, e3, 'm-', label='p=3')
loglog(Nn, f, 'm--', linewidth=1.0)
text(5,1.0e-5, r"$\Delta x^6$", color='m')

gca().set_xlim(0,40)
xticks([2,4,8,16,32], [2,4,8,16,32])
grid()
legend(loc='best')
xlabel('Number of cells')
ylabel(r'$L_2$ Error')

savefig('conv-rate-recovery-dd.png', dpi=150)
show()

