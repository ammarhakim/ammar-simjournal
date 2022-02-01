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

print "RK2, Fixed dx"
dat = pylab.loadtxt("dt-convergence-rk2.dat")
convergenceDt(dat[:,0], dat[:,1])

print "RK3, Fixed dx"
dat = pylab.loadtxt("dt-convergence-rk3.dat")
convergenceDt(dat[:,0], dat[:,1])

print "DG2 convergence"
dat = pylab.loadtxt("dx-convergence-dg-2.dat")
convergence(dat[:,0], dat[:,1])

print "DG3 convergence"
dat = pylab.loadtxt("dx-convergence-dg-3.dat")
convergence(dat[:,0], dat[:,1])

print "DG2 Enstrophy convergence"
dat = pylab.loadtxt("enstrophy-convergence-dg-2.dat")
convergence(dat[:,0], dat[:,1])
