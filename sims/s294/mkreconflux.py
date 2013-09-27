import numpy
from pylab import *
import tables

rc('text', usetex=True)

Lx = 50.0
Ly = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
elcMass = ionMass/25
ionCycl = ionCharge*B0/ionMass

start = 0
end  = 60
nFrame = end-start+1
tm = zeros((nFrame,), float)
flx = zeros((nFrame,), float)

count = 0
for i in range(start, end+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("s294-harris-tenmom_q_%d.h5" % i)
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]
    YI = ny/4

    X = linspace(0, Lx, nx)
    Y = linspace(0, Ly, ny)

    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    
    tm[count] = fh.root.timeData._v_attrs['vsTime']
    flx[count] = dx*sum(abs(q[0:nx,YI,24]))

    count = count+1

def calcDeriv(T, func):
    nt = T.shape[0]-1
    tm = numpy.zeros((nt,), numpy.float)
    vx = numpy.zeros((nt,), numpy.float)

    for i in range(nt):
        tm[i] = 0.5*(T[i+1]+T[i])
        vx[i] = (func[i+1]-func[i])/(T[i+1]-T[i])

    return tm, vx

tmDiff, flxDiff = calcDeriv(tm, flx)

figure(1)
plot(ionCycl*tm, flx, '-ko', label='$\psi$')
legend(loc='best')
title('$d\psi/dt$')
xlabel('Time')
ylabel('$\psi$')

figure(2)
plot(ionCycl*tmDiff, flxDiff, '-ko')
legend(loc='best')
title('$d\psi/dt$')
xlabel('Time')
ylabel('$d\psi/dt$')

show()
