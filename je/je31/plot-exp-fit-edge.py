from pylab import *
import scipy.optimize

style.use('postgkyl.mplstyle')

cfl = 0.1
R = linspace(-3, 3, 100)

# exact exponential fit
def func(an, f0, f1):
    a0 = an[0]; a1 = an[1]
    rhs0 = (exp(a1+a0) - exp(a0-a1))/a1
    rhs1 = ((a1-1)*exp(a0+a1) + (a1+1)*exp(a0-a1))/a1**2

    return rhs0-2*f0, rhs1-2.0/3.0*f1

R1 = linspace(-2.99, 2.99, 100)
gR = R1*0.0
gL = R1*0.0
for i in range(R1.shape[0]):
    aout = scipy.optimize.fsolve(func, [1.0, 0.01], args=(1.0, R1[i]))
    gR[i] = exp(aout[0]+aout[1])
    gL[i] = exp(aout[0]-aout[1])

fig, ax1 = subplots()
xlabel('r = f_1/f_0')

plot(R1, gR, '-r')
plot(R1, gL, '-b')
plot(R1, 1+R1, '--r')
plot(R1, 1-R1, '--b')
plot([R[0],R[-1]], [1/cfl, 1/cfl], 'g--', linewidth=1.0)
xlabel('r')
title('gR/f0 (red) gL/f0 (blue)')



gca().set_ylim([-2,2/cfl])
grid()

savefig('exp-fit-edge.png', dpi=150)
show()

