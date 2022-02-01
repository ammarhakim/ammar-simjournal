#!/usr/bin/env python

import math
import numpy as np
import scipy.optimize as opt
import scipy.special as spc
import scipy.fftpack as fft

def Z(zeta):
    return 1j * np.sqrt(np.pi) * np.exp(-zeta**2) * (1+spc.erf(1j*zeta))
def derZ(zeta):
    return -2*(1 + zeta*Z(zeta))
def der2Z(zeta):
    return -2*(Z(zeta) + zeta*derZ(zeta))

def WeibelDispRel(omega, k, u, vth, c=1, omega_pe=1):
    zeta = omega/k/np.sqrt(2*vth**2)
    return 1.0 - omega_pe**2/(c**2*k**2)*(zeta*Z(zeta)*(1 + u**2/vth**2) + u**2/vth**2) - 2*vth**2/c**2 * zeta**2
def derWeibelDispRel(omega, k, u, vth, c=1, omega_pe=1):
    zeta = omega/k/np.sqrt(2*vth**2)
    return - omega_pe**2/(c**2*k**2) * (1 + u**2/vth**2) * (Z(zeta) + zeta*derZ(zeta))/(k*vth*np.sqrt(2)) - 4*zeta*vth/(c**2*k*np.sqrt(2))

def WeibelDispRelT(omega, k, vthx, vthy, c=1, omega_pe=1):
    zeta = omega/k/np.sqrt(2*vthx**2)
    return 1.0 - omega_pe**2/(c**2*k**2)*(vthy**2/vthx**2*(1+zeta*Z(zeta))-1) - 2*vthx**2/c**2*zeta**2

vt = math.sqrt(0.1)
kv = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

for i in range(len(kv)):
    k = kv[i]
    theor = opt.newton(WeibelDispRel, 0.04j, fprime=derWeibelDispRel, args=(k, 1.5*vt, vt), maxiter=10000)
    print(k, theor.imag)
