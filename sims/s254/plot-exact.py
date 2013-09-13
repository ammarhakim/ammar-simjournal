from pylab import *
from scipy import special

def exactSol(T):
    Nt = T*0.0
    for i in range(T.shape[0]):
        t = T[i]
        t1 = pi*special.erf(1/sqrt(2)*min(6,2*pi/t))
        t2 = t/sqrt(2*pi)*(1-exp(-min(6,2*pi/t)**2/2))
        Nt[i] = t1-t2
    return Nt

def exactEnergySol(T):
    Et = T*0.0
    for i in range(T.shape[0]):
        t = T[i]
        t1 = pi*special.erf(min(6,sqrt(2)*pi/t))
        t2 = t*sqrt(2/pi)*(1-exp((-min(6,2*pi/t)**2/2)))
        Et[i] = t1-t2
    return Et

T = linspace(0.01, 20, 20)

figure(1)
Nt = exactSol(T)
tPtcl = loadtxt("s254-free-stream-bounded_totalPtcl.txt")

plot(T, 2*Nt, 'ro', label='Exact')
plot(tPtcl[:,0], tPtcl[:,1], '-k', label='Gkeyll')
legend()
title('Total Particles')
ylabel('Number of Particles')
xlabel('Time [s]')
savefig('s254-totalPtcl-fs-bounded_tPtcl.png')

figure(2)
Et = exactEnergySol(T)
tPtclEnergy = loadtxt("s254-free-stream-bounded_totalPtclEnergy.txt")

plot(T, 2*Et, 'ro', label='Exact')
plot(tPtclEnergy[:,0], tPtclEnergy[:,1], '-k', label='Gkeyll')
legend()
title('Particle Kinetic Energy')
ylabel('Kinetic Energy')
xlabel('Time [s]')
savefig('s254-totalPtcl-fs-bounded_tEnergy.png')
show()
