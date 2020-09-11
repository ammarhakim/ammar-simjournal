from pylab import *
import postgkyl as pg

style.use('../postgkyl.mplstyle')

komega = pg.GData("iso-weibel-1_frequencies.bp").getValues()
k = komega[:,0]

numk = k.shape[0]
ntot = komega.shape[1]
numEig = int((ntot-3)/2)

omega_r = zeros((numk, numEig), float)
omega_i = zeros((numk, numEig), float)
for i in range(numk):
    omega_r[i,:] = komega[i,3::2]
    omega_i[i,:] = komega[i,4::2]

figure(1)
for i in range(numEig):
    plot(0.1*k, abs(omega_i[:,i]), '.', color='r', markersize=2)
grid()

u = 0.1
# plot the exact solution
wexact = -(sin(arctan2(0,(-sqrt(8*k**2*u**2+k**4+4*k**2+4))+k**2+2)/2.0+0)*sqrt(abs(sqrt(8*k**2*u**2+k**4+4*k**2+4)-k**2-2)))/sqrt(2)

#plot(0.1*k, -wexact, '-r')

xlabel(r'$k u_d/\omega_{pe}$')
ylabel(r'$\gamma/\omega_{pe}$')

savefig('iso-weibel-growth.png', dpi=200)

show()
