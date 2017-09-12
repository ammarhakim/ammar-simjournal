from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

# total change

maf = loadtxt("m4-2d-adv-dg_deltaChange.txt")
saf = loadtxt("../s4/s4-2d-adv-dg_deltaChange.txt")

figure(1)
subplot(2,1,1)
plot(maf[:,0], maf[:,1], label='AL')
plot(saf[:,0], saf[:,1], label='NL')
ylabel(r'$\Delta f$')
legend(loc='best')
grid()

maf = loadtxt("m4-2d-adv-dg_rescaledCells.txt")
saf = loadtxt("../s4/s4-2d-adv-dg_rescaledCells.txt")

subplot(2,1,2)
plot(maf[:,0], maf[:,1], label='AL')
plot(saf[:,0], saf[:,1], label='NL')
xlabel('Time [s]')
ylabel(r'# Cells')
legend(loc='best')
grid()

savefig('s4-m4-df-nc-cmp.png', dpi=150)

show()
