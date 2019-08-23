from pylab import *

cfl = 0.1

def getAlpha(r):
    if r < 2.2:
        return (1+r/3.0)*exp(2.0*r/3.0)
    else:
        return min(1/cfl, 6/(3-r))

def getDgAlpha(r):
    return 1+r
    
def getRn(r):
    al = getAlpha(r)
    return (r-3*cfl*al+6*cfl)/(1-cfl*al)
    
r = linspace(0, 2.39, 100)
rn = r*0.0
alr = r*0.0

for i in range(r.shape[0]):
    alr[i] = getAlpha(r[i])

#figure(1)
#plot(r, alr)
#grid()    
    
for i in range(r.shape[0]):
    rn[i] = getRn(r[i])

figure(2)
plot(r, rn)
grid()
    
show()

