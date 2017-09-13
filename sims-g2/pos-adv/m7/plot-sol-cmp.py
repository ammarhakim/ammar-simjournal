from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

# total change

maf7 = loadtxt("m7-2d-adv-dg_deltaChange.txt")
maf6 = loadtxt("../m6/m6-2d-adv-dg_deltaChange.txt")

mad = loadtxt("m7-2d-adv-dg_density.txt")
N0 = mad[0,1] # initial density

figure(1)
plot(maf7[:,0], maf7[:,1]/N0, label='dt/2')
plot(maf6[:,0], maf6[:,1]/N0, label='dt')
ylabel(r'$\frac{1}{N_0}\Delta f$')
legend(loc='best')
grid()

savefig('m7-m6-df-cmp.png', dpi=150)

e6 = sum(maf6[:,1])
e7 = 0.5*sum(maf7[:,1])
print("Change in error is %g" % (e6/e7))

show()
