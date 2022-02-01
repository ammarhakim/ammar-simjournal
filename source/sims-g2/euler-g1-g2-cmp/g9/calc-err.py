from pylab import *
import postgkyl

u = 1
v = -0.5

def exactRho(X, Y, tm):
    return 1 + 0.2*sin(pi*((X+Y)-tm*(u+v)))
    

d = postgkyl.GData("g9-euler-ds-2d_fluid_1.bp")

q = d.getValues()
tm = d.time
g = d.getGrid()
dx = g[0][1]-g[0][0]
dy = g[1][1]-g[1][0]

X = linspace(g[0][0]+0.5*dx, g[0][-1]-0.5*dx, g[0].shape[0]-1)
Y = linspace(g[1][0]+0.5*dx, g[1][-1]-0.5*dx, g[1].shape[0]-1)
XX, YY = meshgrid(X, Y)

qEx = exactRho(XX, YY, tm)

err = sqrt(dx*dy*sum((q[:,:,0]-qEx)**2))
print(dx, err)

