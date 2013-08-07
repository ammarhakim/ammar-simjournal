from pylab import *
from scipy import special
erf = special.erf

def exactSol(T):
    Nt = T*0.0
    for i in range(T.shape[0]):
        t = T[i]
        # auto-generatef from Maxima
        Nt[i] = pi*erf(1/sqrt(2))+((t+sqrt(2)*pi**(3.0/2.0)*exp(min(6,2*pi/t)**2/2.0)*erf(min(6,2*pi/t)/sqrt(2)))*exp(-min(6,2*pi/t)**2/2.0)-exp((-1.0)/2.0)*(t+sqrt(2)*exp(1.0/2.0)*pi**(3.0/2.0)*erf(1/sqrt(2))))/(sqrt(2)*sqrt(pi))
    return Nt

T = linspace(0.01, 2.0*pi, 20)
Nt = exactSol(T)
tPtcl = loadtxt("s269-free-stream-bounded_totalPtcl.txt")

plot(T, 2*Nt, 'ro', label='Exact')
plot(tPtcl[:,0], tPtcl[:,1], '-k', label='Gkeyll')
legend()
title('Total Particles')
ylabel('Number of Particles')
xlabel('Time [s]')
savefig('s269-totalPtcl-fs-bounded.png')
show()
