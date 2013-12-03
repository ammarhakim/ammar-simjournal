from pylab import *
import tables

dat_ldg = loadtxt("ldg_l2Norm.txt")
dat_rdg = loadtxt("rdg_l2Norm.txt")

solL2 = 3.1399614 # this is L2 norm of initial condition

semilogy(dat_ldg[:,0], sqrt(dat_ldg[:,1]/(dat_ldg[:,0]*solL2)), 'r-', label='LDG')
semilogy(dat_rdg[:,0], sqrt(dat_rdg[:,1]/(dat_rdg[:,0]*solL2)), 'k-', label='RDG')
legend(loc='best')

show()


