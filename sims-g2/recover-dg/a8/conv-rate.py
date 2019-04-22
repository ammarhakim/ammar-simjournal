from pylab import *
import math

style.use('../code/postgkyl.mplstyle')

dat = loadtxt("error-dx.txt")
dx = dat[:,0]
err = dat[:,1]

for i in range(1,dx.shape[0]):
    dxOrder = math.log(err[i-1]/err[i])/log(dx[i-1]/dx[i])
    print("%g %g" % (dx[i], dxOrder))

dx1 = linspace(1.5, 6.0, 10)
err4 = (1/dx1)**4 
    
loglog(1/dx, err)
loglog(dx1, err4, 'r--')
xticks([])
text(3, 0.02, r'$\Delta x^4$')
xlabel(r'$1/\Delta x$')
ylabel('Error')
grid()
show()

    
