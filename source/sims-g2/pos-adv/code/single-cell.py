from pylab import *

style.use('../code/postgkyl.mplstyle')

cfl = 0.1
nSteps = 40

def alpha(f0, f1, x, CFL):
    r = f1/(f0 + 1e-16)
    val = 0.0
    if x > 0:
        if r<2.2:
	    val = max(0, min(1.0/CFL, exp(2*r*x/3)*(1+r*x/3)))
        else:
	    val = min(1.0/CFL, 6/(3-min(3.0, abs(r))))
    else:
        if r>-2.2:
	    val = max(0, min(1.0/CFL, exp(2*r*x/3)*(1+r*x/3)))
        else:
	    val = min(1.0/CFL, 6/(3+min(3.0, abs(r))))
    return val

T = linspace(0, cfl*nSteps, nSteps+1)
f0 = zeros((nSteps+1,), float)
f1 = zeros((nSteps+1,), float)

f0[0] = 1.0
f1[0] = 0.0

def marchCell(f0, f1):
    for i in range(nSteps):
        f0c = f0[i]
        f1c = f1[i]
        alp = alpha(f0c,f1c,1,cfl)
        
        f0[i+1] = f0c*(1-cfl*alp)
        f1[i+1] = f1c - 3*cfl*(alp-2)*f0c

marchCell(f0, f1)

plot(T, f0, label='f0')
plot(T, f1, label='f1')
plot(T, f1/f0, label='f1/f0')
xlabel('Time')
title('Time-history of solution in single cell')
#plot([T[0], T[-1]], [sqrt(2),sqrt(2)],'--')
#gca().set_ylim([0,3])
legend(loc='best')
grid()
savefig('single-cell-hist.png')
show()
