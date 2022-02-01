from pylab import *
import postgkyl as pg

style.use('../postgkyl.mplstyle')

komega = pg.GData("10m-weibel-3_frequencies.bp").getValues()
k = komega[:,0]

numk = k.shape[0]
ntot = komega.shape[1]
numEig = int((ntot-3)/2)

omega_r = zeros((numk, numEig), float)
omega_i = zeros((numk, numEig), float)
for i in range(numk):
    omega_r[i,:] = komega[i,3::2]
    omega_i[i,:] = komega[i,4::2]

lD = sqrt(0.1)
    
figure(1)
for i in range(numEig):
    plot(lD*k, omega_i[:,i], '.', color='k', markersize=2) 

kinDat = loadtxt("kin-growth")
plot(lD*kinDat[:,0], kinDat[:,1], 'o', color='r')
    
#gca().set_xlim([0.0, 0.1])
#gca().set_ylim([0.0, 0.08])
grid()

xlabel(r'$k \lambda_D$')
ylabel(r'$\gamma/\omega_{pe}$')

savefig('10m-weibel-growth.png', dpi=200)

show()
