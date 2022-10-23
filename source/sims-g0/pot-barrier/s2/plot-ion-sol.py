from pylab import *
import postgkyl as pg

m = 1.0
vth = 0.5
deltaPhi = -0.4
charge = 1.0
udrift = 0.6

VR = linspace(-3.0, 3.0, 1024)
vcut = sqrt(-2*charge*deltaPhi/m)
VRc = linspace(vcut, 3.0, 500)
VLc = sqrt( abs(2.0/m*(0.5*m*VRc**2 + charge*deltaPhi)) )

def maxwellian(V, n, u, vth):
    return n/sqrt(2*pi*vth**2)*exp(-(V-u)**2/(2*vth**2))

data = pg.Data("s2-sheath-pot-ion_20.gkyl")
ion_interp = pg.data.GInterpModal(data, 2, "ms")
ion_interp.interpolate(overwrite=True)
pg.data.select(data, z0=0.9, overwrite=True)

# plot Gkyl data
pg.output.plot(data)

plot(VR, maxwellian(VR, 1.0, udrift, vth), "-r")
plot(VRc, maxwellian(VLc, 1.0, udrift, vth), "-k")
xlim([-2.0, 3.0])
grid("on")

show()
