from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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

#majorLocator   = MultipleLocator(0.1)
#majorFormatter = FormatStrFormatter('%g')

def dampOrder1(K):
    return -math.sqrt(math.pi/8)*(1/K**3)*pylab.exp(-0.5/K**2-1.5)

def dampOrder3(K):
    return -math.sqrt(math.pi/8)*(1/K**3-6*K)*pylab.exp(-0.5/K**2-1.5-3*K**2-12*K**4)

dat = loadtxt('ld-damping-rates-elc-osc.txt')
K = pylab.linspace(math.sqrt(0.04), 0.45, 100)
f = figure(2)
sp = subplot(111)
#loglog(K, -dampOrder3(K), '-r', label='O3')
#loglog(K, -dampOrder1(K), '-b', label='O1')
#loglog(dat[:,0], dat[:,1], 'ko')
#semilogy(K, -dampOrder3(K), '-r', label='O3')
#semilogy(K, -dampOrder1(K), '-b', label='O1')
semilogy(dat[:,0], dat[:,1], 'ko')
datEx = loadtxt('exact-ld-damping-rates-elc-osc.txt')
semilogy(datEx[:,0], -datEx[:,1], 'm-', label='Exact')
pylab.legend(loc='upper left')
#sp.xaxis.set_major_locator(majorLocator)
#sp.xaxis.set_major_formatter(majorFormatter)
xlabel(r'$k\lambda_D$')
ylabel('Damping rate')
pylab.axis('tight')
savefig('ld-damping-rates-elc-osc.png')
close()

exactRates = numpy.array([.007780445578799889, .01844791241964087, .03240172691301003, 
                          .1533594669096014, 0.461918799798943, .8513304586905912, 3.966650346524154])
print "Percent errors", -100*(dat[:,1]-exactRates)/exactRates
