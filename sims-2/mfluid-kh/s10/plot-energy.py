import matplotlib
from pylab import *
import tables
import math
import tftenmoment

totEm = loadtxt("s10-10m-karim-kh_totEm.txt")
totBxBy = loadtxt("s10-10m-karim-kh_totBxBy.txt")
totIe_e = loadtxt("s10-10m-karim-kh_totIe_e.txt")
totKe_e = loadtxt("s10-10m-karim-kh_totKe_e.txt")
totIe_i = loadtxt("s10-10m-karim-kh_totIe_i.txt")
totKe_i = loadtxt("s10-10m-karim-kh_totKe_i.txt")

T = linspace(0, 500, totEm.shape[0])

figure(1)

subplot(3,2,1)
plot(T, totEm, 'r-')
title('EM')
subplot(3,2,2)
totE = totEm+totKe_e+totKe_i+totIe_e+totIe_i
ylabel('Total')
plot(T, totE/totE[0], 'r-')

subplot(3,2,3)
plot(T, totKe_e/totKe_e[0], 'r-')
ylabel('KE_e')
subplot(3,2,4)
plot(T, totIe_e/totIe_e[0], 'r-')
ylabel('IE_e')

subplot(3,2,5)
plot(T, totKe_i/totKe_i[0], 'r-')
ylabel('KE_i')
subplot(3,2,6)
plot(T, totIe_i/totIe_i[0], 'r-')
ylabel('IE_i')

print("Drop in total energy", (totE[-1]-totE[0])/totE[0]*100)
savefig('s10-energy-profs.png')

def relDiff(v):
    return (v-v[0])/totE[0]

figure(2)
plot(T, relDiff(totIe_e), '-b', label='dE(P_e)')
plot(T, relDiff(totIe_i), '-k', label='dE(P_i)')
plot(T, relDiff(totBxBy), '-r', label='dE(B)')
plot(T, relDiff(totKe_i), '-m', label='dE(U_i)')
legend(loc='best')
savefig('s10-rel-energy-profs.png')

show()
