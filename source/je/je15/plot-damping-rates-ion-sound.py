from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import numpy
import math

figure(1)
rc('text', usetex=True)
dat = loadtxt('damping-rates-ion-sound.txt')
semilogy(dat[:,0], dat[:,1]/(math.sqrt(2)*1.0*0.5), 'ko')
xlabel(r'$T = T_i/T_e$')

datEx = loadtxt('exact-ld-damping-rates-ion-acoustic.txt')
semilogy(datEx[:,0], -datEx[:,1], 'm-', label='Exact')

ylabel('Normalized Damping Rate ($\gamma/\sqrt{2} v_t k$)')
pylab.axis('tight')
savefig('damping-rates-ion-sound.png')
show()
close()
