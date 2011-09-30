import pylab
import tables
import math

dat = pylab.loadtxt("fdtd-convergence.dat")
dx = dat[:,0]
err75 = dat[:,1]
err150 = dat[:,2]

fig = pylab.figure(1)
pylab.loglog(dx, err75, 'r-')
pylab.loglog(dx, err75, 'ro')
pylab.loglog(dx, err150, 'k-')
pylab.loglog(dx, err150, 'ko')
pylab.xlabel("Grid spacing")
pylab.ylabel("Average error")
pylab.axis('tight')
pylab.title('Error for wave-propagation scheme')
pylab.show()
pylab.savefig("tm-wave-convergence.png")

print "Order of convergence t=75"
for i in range(1, dx.shape[0]):
    od = (pylab.log(err75[i]) - pylab.log(err75[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
    print od

print "Order of convergence t=150"
for i in range(1, dx.shape[0]):
    od = (pylab.log(err150[i]) - pylab.log(err150[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
    print od
