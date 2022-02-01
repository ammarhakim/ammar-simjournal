from pylab import *

def loadData(sid):
    tm = loadtxt("../s%d/s%d-recon-tmPoints" % (sid, sid))
    flux = loadtxt("../s%d/s%d-recon-reconFlux" % (sid, sid))
    return tm, flux

def calcGrad(x, fx):
    xmid = zeros((x.shape[0]-1,), float)
    dfx = zeros((x.shape[0]-1,), float)
    for i in range(x.shape[0]-1):
        xmid[i] = 0.5*(x[i+1]+x[i])
        dfx[i] = (fx[i+1]-fx[i])/(x[i+1]-x[i])

    return xmid, dfx

# load data
tm5, f5 = loadData(433)
tm10, f10 = loadData(435)
tm15, f15 = loadData(437)
tm25, f25 = loadData(436)

# plot on same figure
figure(1)
plot(tm5, f5, 'r-', label='5di')
plot(tm10, f10*(5.0/10.0), 'k-', label='10di')
plot(tm15, f15*(5.0/15.0), 'b-', label='15di')
plot(tm25, f25*(5.0/25.0), 'g-', label='25di')
legend(loc='upper left')
title('Normalized reconnected flux')

# compute rates
tm5, r5 = calcGrad(tm5, f5)
tm10, r10 = calcGrad(tm10, f10*(5.0/10.0))
tm15, r15 = calcGrad(tm15, f15*(5.0/15.0))
tm25, r25 = calcGrad(tm25, f25*(5.0/25.0))

# plot on same figure
figure(2)
plot(tm5, r5, 'r-', label='5di')
plot(tm10, r10, 'k-', label='10di')
plot(tm15, r15, 'b-', label='15di')
plot(tm25, r25, 'g-', label='25di')
legend(loc='upper left')
title('Reconnection rate')

show()
