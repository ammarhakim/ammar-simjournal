from pylab import *
import pylab
import numpy
import math

figure(1)
rc('text', usetex=True)
dat = loadtxt('damping-rates-elc-osc.txt')
semilogy(0.5*sqrt(dat[:,0]), dat[:,1], '-ro')
xlabel(r'$k\lambda_D$')
ylabel('Damping rate')
savefig('damping-rates-elc-osc.png')
close()

def dampOrder1(K):
    return -math.sqrt(math.pi/8)*(1/K**3)*pylab.exp(-0.5/K**2-1.5)

def dampOrder3(K):
    return -math.sqrt(math.pi/8)*(1/K**3-6*K)*pylab.exp(-0.5/K**2-1.5-3*K**2-12*K**4)

dat = loadtxt('ld-damping-rates-elc-osc.txt')
K = pylab.linspace(math.sqrt(0.04), math.sqrt(0.12), 100)
figure(2)
semilogy(K**2, -dampOrder3(K), '-r', label='O3')
semilogy(K**2, -dampOrder1(K), '-b', label='O1')
semilogy(dat[:,0]**2, dat[:,1], 'ko')
pylab.legend(loc='upper left')
xlabel(r'$(k\lambda_D)^2$')
ylabel('Damping rate')
savefig('ld-damping-rates-elc-osc.png')
close()

print "Percent errors", -100*(dat[:,1]+dampOrder3(dat[:,0]))/dampOrder3(dat[:,0])
