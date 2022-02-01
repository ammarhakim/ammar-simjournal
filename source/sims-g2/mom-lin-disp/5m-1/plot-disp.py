from pylab import *
import postgkyl as pg

style.use('../postgkyl.mplstyle')

def getReal(fn):
    komega = pg.GData(fn).getValues()
    k = komega[:,0]

    numk = k.shape[0]
    ntot = komega.shape[1]
    numEig = int((ntot-3)/2)

    # collect all real-part of eigenvalues
    omega_r = zeros((numk, numEig), float)
    for i in range(numk):
        omega_r[i,:] = komega[i,3::2]
        omega_r.sort()

    return k, omega_r

k_5m, omega_r_5m = getReal("5m-1-waves_frequencies.bp")
k_iso, omega_r_iso = getReal("../iso-1/iso-1-waves_frequencies.bp")

figure(1)
for i in range(16):
    plot(k_5m, abs(omega_r_5m[:,i]), '.', color='r', markersize=2)

for i in range(14):
    plot(k_iso, abs(omega_r_iso[:,i]), '.', color='k', markersize=2)

grid()
xlabel("k")
ylabel("$\omega_r$")
title("Cold (black) and five-moment (red) dispersion relations")
savefig("iso-5m-cmp-cold.png", dpi=200, bbox='tight')
    
show()

