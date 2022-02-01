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

k_5m, omega_r_5m = getReal("../5m-1/5m-1-waves_frequencies.bp")
k_10m, omega_r_10m = getReal("../10m-1/10m-1-waves_frequencies.bp")

figure(1)
for i in range(16):
    plot(k_5m, abs(omega_r_5m[:,i]), '.', color='r', markersize=2)

for i in range(14):
    plot(k_10m, abs(omega_r_10m[:,i]), '.', color='k', markersize=2)

grid()
xlabel("k")
ylabel("$\omega_r$")

gca().set_xlim([0.0, 3.0])
gca().set_ylim([0.0, 1.0])

#text(0.01, 0.05 + 0.01, "$\omega_{ci}$")
#text(0.01, 0.1 + 0.01, "$2 \omega_{ci}$")

title("Ten-moment (black) and five-moment (red)")
#savefig("10m-5m-cmp-ion.png", dpi=200, bbox='tight')
    
show()

