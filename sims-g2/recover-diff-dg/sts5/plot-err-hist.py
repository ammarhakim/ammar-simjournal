from pylab import *

style.use('../postgkyl.mplstyle')

st = loadtxt("sts5-recovery-ss-diff_err.txt")

numStages = int(loadtxt("numStages"))

def calcSlope(s, e, ns, step, dat):
    return (log(dat[s])-log(dat[e]))/(s-e)/ns

nt = st[:,0]
dt = st[:,1]

semilogy(numStages*nt, dt, 'r-')

startIdx = nt.shape[0]//2
endIdx = nt.shape[0]-1
pt = calcSlope(startIdx, endIdx, numStages, nt, dt)
print(pt, pt*numStages)

semilogy(numStages*nt, 1e-1*exp(pt*numStages*nt), 'r--')

grid()
xlabel('RK-Stages (s)')
ylabel(r'$L_2$ Error')
title("Num Stages = %d." % numStages)

savefig('ss5-sts5-error.png', dpi=200)

show()

