from pylab import *
import postgkyl as pg

mass = 1.0
charge = 1.0
density = 1.0
temperature = 0.1
Bx = 1.0
Bz = 0.75

Bmag = math.sqrt(Bx**2 + Bz**2)
omegac = charge*Bmag/mass
vthe = math.sqrt(temperature/mass)
rhoe = vthe/omegac

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
    plot(rhoe*k_5m, abs(omega_r_5m[:,i]), '.', color='r', markersize=2)

for i in range(14):
    plot(rhoe*k_10m, abs(omega_r_10m[:,i]), '.', color='k', markersize=2)

grid()
xlabel(r"$k\rho_e$")
ylabel("$\omega_r$")

gca().set_xlim([rhoe*0.0, rhoe*3.0])
gca().set_ylim([0.5, 3.5])

text(rhoe*0.01, 1.25+0.1, "$\omega_{ce}$")
text(rhoe*0.01, 2.5+0.1, "$2 \omega_{ce}$")

title("Ten-moment (black) and five-moment (red)")
savefig("10m-5m-cmp-elc.png", dpi=200, bbox='tight')
    
show()

