import pylab
import math

def convergence(dx, err):
    for i in range(1, dx.shape[0]):
        od = (pylab.log(err[i]) - pylab.log(err[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
        print od

print "Dirichlet BCs"
dat = pylab.loadtxt("convergence-1d-d.dat")
convergence(dat[:,0], dat[:,1])

print "Neumann BCs"
dat = pylab.loadtxt("convergence-1d-n.dat")
convergence(dat[:,0], dat[:,1])

print "Two dimensions"
dat = pylab.loadtxt("convergence-2d.dat")
convergence(dat[:,0], dat[:,1])

print "One D, Order 3"
dat = pylab.loadtxt("convergence-1d-o3.dat")
convergence(dat[:,0], dat[:,1])

print "One D, Order 4"
dat = pylab.loadtxt("convergence-1d-o4.dat")
convergence(dat[:,0], dat[:,1])

print "Two dimensions, Order 3"
dat = pylab.loadtxt("convergence-2d-o3.dat")
convergence(dat[:,0], dat[:,1])

print "Two dimensions, Periodic"
dat = pylab.loadtxt("convergence-periodic.dat")
convergence(dat[:,0], dat[:,1])
