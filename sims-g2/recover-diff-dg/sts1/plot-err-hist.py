from pylab import *

style.use('../postgkyl.mplstyle')

st = loadtxt("sts1-recovery-ss-diff_err.txt")
ss = loadtxt("../ss1/ss1-recovery-ss-diff_err.txt")

numStages = 10

nt = st[:,0]
dt = st[:,1]

ns = ss[:,0]
ds = ss[:,1]

semilogy(numStages*nt, dt, color='red')
semilogy(3*ns, ds, color='black')
grid()
xlabel('RK-Stages')
ylabel(r'$L_2$ Error')
title("Num Stages = %d. Speedup = %g" % (numStages, (3*ns.shape[0]/(numStages*nt.shape[0]))))

savefig('ss1-sts1-error.png', dpi=200)

show()

