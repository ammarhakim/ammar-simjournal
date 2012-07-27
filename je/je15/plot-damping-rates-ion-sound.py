from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pylab
import numpy
import math

figure(1)
rc('text', usetex=True)
dat = loadtxt('damping-rates-ion-sound.txt')
ndr = dat[:,1]/(0.5*sqrt(1/dat[:,0]))
semilogy(dat[:,0], ndr, '-ro')
xlabel(r'$T = T_i/T_e$')
ylabel('Normalized Damping Rate (\gamma/c_ek)')
savefig('damping-rates-ion-sound.png')
close()
