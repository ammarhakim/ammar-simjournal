import numpy
import pylab
import tables
import optparse
import math

gasGamma = 2.0
rho0 = 1.0
p0 = 1.0
U0 = 1e-8
nModes = 5000
cs0 = math.sqrt(gasGamma*p0/rho0)

def exactSol(X, t, Lambda):
    Bz = 1.0
    wc = Lambda*Bz
    u1 = 0*X
    U0 = 1e-8
    for n in range(nModes+1):
        kn = 2*math.pi*(2*n+1)
        wn = math.sqrt(kn*kn*cs0*cs0 + wc*wc)
        u1 = u1 - U0/(2*n+1)*numpy.sin(kn*X + wn*t)

    return u1

def main():
    dxExact = (1.0-0.0)/10000
    Xex = pylab.linspace(0.5*dxExact, 1-0.5*dxExact, 10000)
    uEx = exactSol(Xex, 1000, 10)
    pylab.plot(Xex, uEx, '-k')
    pylab.xlabel("X")
    pylab.ylabel("X-component of velocity")
    pylab.savefig("s41-sqpulse-exact.png")

if __name__ == '__main__': main()

