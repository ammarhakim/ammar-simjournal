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


print "DG2 Energy convergence"
dat = pylab.loadtxt("total-energy-convergence-ds.dat")
convergence(dat[:,0], dat[:,1])

print "DG2 Energy convergence (CF)"
dat = pylab.loadtxt("total-energy-convergence-ds-cf.dat")
convergence(dat[:,0], dat[:,1])

print "DG2 Enstrophy convergence (CF)"
dat = pylab.loadtxt("total-enstrophy-convergence-ds-cf.dat")
convergence(dat[:,0], dat[:,1])
