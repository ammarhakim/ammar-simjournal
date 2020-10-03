from pylab import *

style.use('../postgkyl.mplstyle')

st = loadtxt("sts3-recovery-ss-diff_err.txt")
ss = loadtxt("../ss3/ss3-recovery-ss-diff_err.txt")

numStages = int(loadtxt("numStages"))

def calcSlope(s, e, ns, step, dat):
    return (log(dat[s])-log(dat[e]))/(s-e)/ns

nt = st[:,0]
dt = st[:,1]

ns = ss[:,0]
ds = ss[:,1]

semilogy(numStages*nt, dt, 'r-')
#semilogy(3*ns, ds, 'k-')

startIdx = nt.shape[0]//2
endIdx = nt.shape[0]-1
pt = calcSlope(startIdx, endIdx, numStages, nt, dt)
ps = calcSlope(1500, 2500, 3, ns, ds)

print(pt, ps)

semilogy(numStages*nt, 1e-1*exp(pt*numStages*nt), 'r--')
#text(600, 1e-3, '$e^{-0.015 s}$')

#semilogy(3*ns, 1.5e-4*exp(ps*3*ns), 'k--')
#text(4000, 5e-6, '$e^{-0.001 s}$')

grid()
xlabel('RK-Stages (s)')
ylabel(r'$L_2$ Error')
title("Num Stages = %d. Speedup = %g" % (numStages, (3*ns.shape[0]/(numStages*nt.shape[0]))))

savefig('ss3-sts3-error.png', dpi=200)

show()

