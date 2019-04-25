from pylab import *
import math

style.use('../code/postgkyl.mplstyle')

dat = loadtxt("error-dx.txt")
N = dat[:,0]*1.0
err = dat[:,1]
err1 = err- err[-1] # this gets rid of dt errors assuming dx is very small

for i in range(1,N.shape[0]-1):
    dxOrder = math.log(err1[i-1]/err1[i])/log(N[i]/N[i-1])
    print("%g %g" % (1/N[i], dxOrder))

fig, ax = plt.subplots(1,1)
    
ax.loglog(N[:-1], err1[:-1])
ax.set_xticks([2,4,8,16])
ax.set_xlim(0,20)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel(r'$N$')
ax.set_ylabel('Error')
grid()
show()
