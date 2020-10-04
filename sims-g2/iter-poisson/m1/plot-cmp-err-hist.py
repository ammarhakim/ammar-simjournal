from pylab import *

style.use('../postgkyl.mplstyle')

fno = loadtxt("m1-iter-poisson_err.txt")
fit = loadtxt("m1-iter-poisson-extra_err.txt")

gno = loadtxt("../s1/s1-iter-poisson_err.txt")
git = loadtxt("../s1/s1-iter-poisson-extra_err.txt")

numStages_f = int(loadtxt("m1-iter-poisson_numStages"))
numStages_g = int(loadtxt("../s1/s1-iter-poisson_numStages"))

def plotSlope(fno, fit, gno, git):
    semilogy(numStages_f*fno[:,0], fno[:,1], 'r-')
    semilogy(numStages_f*fit[:,0], fit[:,1], 'k-')

    semilogy(numStages_g*gno[:,0], gno[:,1], 'r--')
    semilogy(numStages_g*git[:,0], git[:,1], 'k--')

    grid()
    xlabel('RK-Stages (s)')
    ylabel(r'$L_2$ Error')
    title("Num Stages = (%d, %d). Dashed: RKL2. Solid: RKL1" % (numStages_f, numStages_g))
    
    savefig('s1-m1-error.png', dpi=200)

    show()

plotSlope(fno, fit, gno, git)
