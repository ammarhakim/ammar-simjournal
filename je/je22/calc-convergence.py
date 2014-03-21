import pylab
import math

def convergence(dx, err):
    for i in range(1, dx.shape[0]):
        od = (pylab.log(err[i]) - pylab.log(err[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
        print od

print "Dimensional splitting"
dat = pylab.loadtxt("conv-smooth-euler.txt")
convergence(dat[:,0], dat[:,1])

