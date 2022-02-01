from pylab import *
import postgkyl as pg

style.use('../postgkyl.mplstyle')

komega = pg.GData("iso-4-buneman_frequencies.bp").getValues()
k = komega[:,0]

numk = k.shape[0]
ntot = komega.shape[1]
numEig = int((ntot-3)/2)

omega_r = zeros((numk, numEig), float)
omega_i = zeros((numk, numEig), float)
for i in range(numk):
    omega_r[i,:] = komega[i,3::2]
    omega_i[i,:] = komega[i,4::2]

omega_growing_r = zeros((numk, numEig), float)
for i in range(numk):
    for n in range(numEig):
        if omega_i[i,n] < 1e-3:
            omega_r[i,n] = -1 # zap these so the are not seen

figure(1)
for i in range(numEig):
    plot(k, omega_r[:,i], '.', color='r', markersize=2)

gca().set_ylim([0,2.5])
grid()

show()
