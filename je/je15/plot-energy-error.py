from pylab import *
import math

dat2 = loadtxt("momentum-error-dx.txt")
dat3 = loadtxt("momentum-error-o3-dx.txt")

loglog(dat2[:,0], dat2[:,1], 'r-o')
loglog(dat3[:,0], dat3[:,1], 'k-o')
axis('tight')
ylabel('Error in Energy')
xlabel('Time-Step')
#show()
savefig('dg-o2-o3-energy-conservation-errors.png')
