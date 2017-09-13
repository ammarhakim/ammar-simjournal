from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

# total change

maf = loadtxt("m6-2d-adv-dg_deltaChange.txt")
saf = loadtxt("../s6/s6-2d-adv-dg_deltaChange.txt")

mad = loadtxt("m6-2d-adv-dg_density.txt")
N0 = mad[0,1] # initial density

figure(1)
subplot(2,1,1)
plot(maf[:,0], maf[:,1]/N0, label='AL')
plot(saf[:,0], saf[:,1]/N0, label='NL')
ylabel(r'$\frac{1}{N_0}\Delta f$')
legend(loc='best')
grid()

maf = loadtxt("m6-2d-adv-dg_rescaledCells.txt")
saf = loadtxt("../s6/s6-2d-adv-dg_rescaledCells.txt")

subplot(2,1,2)
plot(maf[:,0], maf[:,1], label='AL')
plot(saf[:,0], saf[:,1], label='NL')
xlabel('Time [s]')
ylabel(r'# Cells')
legend(loc='best')
grid()

savefig('s6-m6-df-nc-cmp.png', dpi=150)

show()
