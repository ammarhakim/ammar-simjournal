from pylab import *
import postgkyl as pg

style.use('../postgkyl.mplstyle')

komega = pg.GData("iso-1-waves_frequencies.bp").getValues()
k = komega[:,0]

numk = k.shape[0]
ntot = komega.shape[1]
numEig = int((ntot-3)/2)

# collect all real-part of eigenvalues
omega_r = zeros((numk, numEig), float)
for i in range(numk):
    omega_r[i,:] = komega[i,3::2]
    omega_r.sort()

figure(1)
plot(k, omega_r[:,13], '.', markersize=2)
plot(k, omega_r[:,12], '.', markersize=2)
plot(k, omega_r[:,11], '.', markersize=2)
plot(k, omega_r[:,10], '.', markersize=2)
plot(k, omega_r[:,9], '.', markersize=2)

grid()
xlabel("k")
ylabel("$\omega_r$")
savefig("iso-cold-waves.png", dpi=200, bbox='tight')

show()

