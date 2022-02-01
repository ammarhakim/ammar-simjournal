from numpy import *
from scipy import special
from scipy import optimize
import matplotlib.pyplot as plt
import math
from matplotlib import rc

beta_e_val = 10.0

def plasmaDisp(z):
    return 1j*sqrt(pi)*exp(-z**2)*(1+special.erf(1j*z))

def derivPlasmaDisp(z):
    return -2*(1+z*plasmaDisp(z))

def eps(z,kPerpRho):
    return (2*beta_e_val*z**2-1)*(1+z*plasmaDisp(z)) - kPerpRho**2

def derivEps(z,kPerpRho):
    return (4*beta_e_val*z)*(1+z*plasmaDisp(z)) + (2*beta_e_val*z**2-1)*(plasmaDisp(z)+z*derivPlasmaDisp(z))

# Number of points to calculate damping rate at
nPoints = 1000;
beta_e_list = [10]
kperp_list = linspace(0.01, 1.0, nPoints)
freqList = zeros(nPoints);
dampList = zeros(nPoints);
approxFreqList = zeros(nPoints);

kPerpRho = kperp_list[0]
kPar = 0.5
# Initial guess for z0 = omega/(k*sqrt(2)*vTe) using approximate expression for wave frequency
z0 = kPar/sqrt(beta_e_list[0]+kPerpRho**2)

for index, kp in enumerate(kperp_list):
    z0 = optimize.newton(eps, z0, derivEps, (kp,), 1e-4, 10000)

    freqList[index] = z0.real;
    dampList[index] = z0.imag;    

#plt.figure(1)
#plt.plot(beta_e_list, freqList, 'g-', label='Better Approx')
#plt.figure(2)
#plt.plot(beta_e_list, dampList, 'b-', label='Approx')
#plt.show()

# write to file
fp = open("KAW-beta-10-rates.txt", "w")
for i in range(nPoints):
    fp.writelines("%g %g %g\n" % (kperp_list[i]**2, freqList[i], dampList[i]) )
fp.close()

