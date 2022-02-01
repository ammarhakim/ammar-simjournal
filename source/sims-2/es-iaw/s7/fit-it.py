from pylab import *
import postgkyl.tools
import scipy.optimize

def lowPass(x, dt, C):
    y = 0*x # filtered signal
    alpha = dt/(C+dt)
    y[0] = x[0]
    for i in range(1,x.shape[0]):
        y[i] = alpha*x[i] + (1-alpha)*y[i-1]
    return y

dat = loadtxt("s7-es-iaw_phiInCell.dat")
T = dat[:,0]
E = dat[:,1]
dt = T[1]-T[0]
wH = 50.25 # computed from Maxima

# FFT filter it, removing high-frequency
fE = real(postgkyl.tools.butterFiltering(E/E[0], dt, 0.35*0.9*wH/(2*pi)))
plot(T, fE)
grid()

figure(2)
# pick range for analysis
tLo = T.searchsorted(3.0)
tHi = T.searchsorted(24.0)

T1 = T[tLo:tHi]
fE1 = fE[tLo:tHi]
fE1l = lowPass(fE1, T1[1]-T1[0], 0.15)

plot(T1+0.1, fE1, 'r-')
plot(T1, fE1l, 'k-')
axis('tight')
grid()

def findMax(t1, t2):
    fv = fE1l[T1.searchsorted(t1):T1.searchsorted(t2)]
    am = fv.argmax()
    return T1[T1.searchsorted(t1):T1.searchsorted(t2)][am], fv[am]

t1, v1 = findMax(6.0, 7.5)
t2, v2 = findMax(10.0, 14.0)
t3, v3 = findMax(17.0, 24.0)

# plot these
plot([t1], [v1], 'bo')
plot([t2], [v2], 'bo')
plot([t3], [v3], 'bo')

# compute best fit
def func(an):
    a0 = an[0]
    rhs0 = exp(a0*(t2-t1)) - v2/v1
    rho1 = exp(a0*(t3-t1)) - v3/v1
    return rhs0, rho1

aout = scipy.optimize.fsolve(func, [0.0])
print(aout[0])
plot(T1, v1*exp(aout[0]*(T1-t1)), 'm-')

show()
