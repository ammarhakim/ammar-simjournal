from pylab import *
import postgkyl

style.use('../code/postgkyl.mplstyle')

def plotLines(q, X, cs):
    mu = linspace(-1, 1, 10)
    for i in range(X.shape[0]-1):
        xx = linspace(X[i], X[i+1], 10)
        ff = 0.5*q[i,0] + sqrt(3)/2.0*mu*q[i,1]
        plot(xx, ff, cs)
        plot(0.5*(xx[0]+xx[-1]), 0.5*q[i,0], cs+'o')

        #print(q[i,1]/(sqrt(3)*q[i,0]))

d = postgkyl.GData("m9-2d-adv-dg_distf_0.bp")
q = d.q[:,0,:]
X = linspace(0, 10.0, 10)

d = postgkyl.GData("m9-2d-adv-dg_distf_1.bp")
q1 = d.q[:,0,:]

figure(1)
plotLines(q, X, 'r')
plotLines(q1, X, 'k')
grid()
title('Red: Initial. Black: t=1.0')
savefig('delta-cell.png')
show()


