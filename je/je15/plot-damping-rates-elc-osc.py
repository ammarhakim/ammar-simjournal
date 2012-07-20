from pylab import *

rc('text', usetex=True)
dat = loadtxt('damping-rates-elc-osc.txt')
semilogy(0.5*sqrt(dat[:,0]), dat[:,1], '-ro')
xlabel(r'$k\lambda_D$')
ylabel('Damping rate')
savefig('damping-rates-elc-osc.png')
