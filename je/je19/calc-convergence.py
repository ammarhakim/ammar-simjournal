import pylab
import math

def convergence(dx, err):
    for i in range(1, dx.shape[0]):
        od = (pylab.log(err[i]) - pylab.log(err[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
        print od

print "polyOrder 1 1D"
dat = pylab.loadtxt("dg-1d-p1-err.txt")
convergence(dat[:,0], dat[:,1])

print "polyOrder 2 1D"
dat = pylab.loadtxt("dg-1d-p2-err.txt")
convergence(dat[:,0], dat[:,1])

print "polyOrder 1 2D"
dat = pylab.loadtxt("dg-2d-p1-err.txt")
convergence(dat[:,0], dat[:,1])
