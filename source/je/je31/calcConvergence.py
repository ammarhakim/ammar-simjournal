import pylab
import math

def convergence(dx, err):
    for i in range(1, dx.shape[0]):
        od = (pylab.log(err[i]) - pylab.log(err[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
        print od

print "Naive DG L2"
dat = pylab.loadtxt("s-errors.txt")
convergence(dat[:,0], dat[:,1])

print "AL DG L2"
dat = pylab.loadtxt("m-errors.txt")
convergence(dat[:,0], dat[:,1])

print "Naive DG PROJ"
dat = pylab.loadtxt("s-errors.txt")
convergence(dat[:,0], dat[:,2])

print "AL DG PROJ"
dat = pylab.loadtxt("m-errors.txt")
convergence(dat[:,0], dat[:,2])
