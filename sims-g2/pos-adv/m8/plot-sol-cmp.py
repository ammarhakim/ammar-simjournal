from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

# total change

maf8 = loadtxt("m8-2d-adv-dg_deltaChange.txt")
maf7 = loadtxt("../m7/m7-2d-adv-dg_deltaChange.txt")
maf6 = loadtxt("../m6/m6-2d-adv-dg_deltaChange.txt")

mad = loadtxt("m8-2d-adv-dg_density.txt")
N0 = mad[0,1] # initial density

figure(1)

mm = interp(maf6[:,0], maf8[:,0], maf8[:,1])
plot(maf6[:,0], maf6[:,1]/mm, label=r'$4 \Delta t$')

mm = interp(maf7[:,0], maf8[:,0], maf8[:,1])
plot(maf7[:,0], maf7[:,1]/mm, label=r'$2 \Delta t$')

plot(maf8[:,0], maf8[:,1]/maf8[:,1], label=r'$\Delta t$')

gca().set_xlim([0,0.06])
ylabel(r'$\Delta f/\Delta f_b$')
xlabel('Time [s]')
legend(loc='best')
grid()

savefig('m8-m7-m6-df-cmp.png', dpi=150)

e6 = sum(maf6[:,1])
e7 = 0.5*sum(maf7[:,1])
e8 = 0.25*sum(maf8[:,1])
print("Change in errors are %g and %g" % (e6/e7, e7/e8))

show()
