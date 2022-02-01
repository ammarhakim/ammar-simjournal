from pylab import *

import pylab as plt
plt.style.use('../paper.mplstyle')

emDat = loadtxt("s2-es-iaw_emEnergy.dat")
ionDat = loadtxt("s2-es-iaw_totalPtclEnergyIon.dat")
elcDat = loadtxt("s2-es-iaw_totalPtclEnergyElc.dat")

T = emDat[:,0]

figure(1)
plot(T, emDat[:,1], 'r-')
title('ES Energy')

figure(2)
plot(T, ionDat[:,1], 'r-')
title('Ion particle Energy')

figure(3)
plot(T, elcDat[:,1], 'r-')
title('Electron particle Energy')

figure(4)
eTot = emDat[:,1]+ionDat[:,1]+elcDat[:,1]
plot(T, eTot, 'r-')

show()
