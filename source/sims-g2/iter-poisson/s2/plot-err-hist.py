from pylab import *

style.use('../postgkyl.mplstyle')

fno = loadtxt("s2-iter-poisson_err.txt")
fit = loadtxt("s2-iter-poisson-extra_err.txt")

numStages = int(loadtxt("s2-iter-poisson_numStages"))

def calcSlope(s, e, ns, step, dat):
    return (log(dat[s])-log(dat[e]))/(s-e)/ns

def plotSlope(fno, fit):
    semilogy(numStages*fno[:,0], fno[:,1], 'r-')
    semilogy(numStages*fit[:,0], fit[:,1], 'k-')

    grid()
    xlabel('RK-Stages (s)')
    ylabel(r'$L_2$ Error')
    title("Num Stages = %d" % numStages)
    
    savefig('s2-error.png', dpi=200)

    show()

plotSlope(fno, fit)    
