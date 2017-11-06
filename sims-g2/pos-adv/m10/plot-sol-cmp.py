from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

# lineouts
d = postgkyl.GData("m10-2d-adv-dg_distf_0.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, f0 = dg.interpolate(0)

d = postgkyl.GData("m10-2d-adv-dg_distf_1.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, mf1 = dg.interpolate(0)

d = postgkyl.GData("../s10/s10-2d-adv-dg_distf_1.bp")
dg = postgkyl.GInterpModal(d, 1, 'ms')
X, sf1 = dg.interpolate(0)

figure(1)
ny = f0.shape[1]
plot(X[0], f0[:,ny/2], 'k', label='EX')
plot(X[0], mf1[:,ny/2], label='AL')
plot(X[0], sf1[:,ny/2], label='NL')
xlabel('X')
ylabel('f(X)')
legend(loc='best')
grid()

savefig('s10-m10-cmp.png', dpi=150)

show()
