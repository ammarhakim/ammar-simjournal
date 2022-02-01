import pylab
import math

def convergence(dx, err):
    for i in range(1, dx.shape[0]):
        od = (pylab.log(err[i]) - pylab.log(err[i-1]))/(pylab.log(dx[i]) - pylab.log(dx[i-1]))
        print od

def convergenceDt(dx, err):
    for i in range(2, dx.shape[0]):
        od = pylab.log(err[i-1]/err[i])/pylab.log(dx[i-1]/dx[i])
        print od 

print "RB DG2 convergence"
dat = pylab.loadtxt("rb-dx-convergence-dg-2.dat")
convergence(dat[:,0], dat[:,1])

print "AD DG2 convergence"
dat = pylab.loadtxt("ad-dx-convergence-dg-2.dat")
convergence(dat[:,0], dat[:,1])

print "AD DG2 3P convergence"
dat = pylab.loadtxt("ad-dx-convergence-dg-2-3p.dat")
convergence(dat[:,0], dat[:,1])
