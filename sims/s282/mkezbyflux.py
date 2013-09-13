import numpy
from pylab import *
import tables

rc('text', usetex=True)

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
mu0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*elcMass*n0)
ionCycl = ionCharge*B0/ionMass

start = 10
end  = 40
nFrame = end-start+1
tm = zeros((nFrame,), float)
elcX = zeros((nFrame,), float)
elcO = zeros((nFrame,), float)
flx = zeros((nFrame,), float)
psiCont = zeros((nFrame,), float)

count = 0
for i in range(start, end+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("s282-gem-tenmom_q_%d.h5" % i)
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]

    X = linspace(-LX/2, LX/2, nx)
    Y = linspace(-LY/2, LY/2, ny)

    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    
    psiB_Up = q[:,ny/2+1,27]
    psiB_Dn = q[:,ny/2-1,27]
    psiDiff = (psiB_Up-psiB_Dn)/(2*dy)

    tm[count] = fh.root.timeData._v_attrs['vsTime']
    flx[count] = dx*sum(q[0:nx/2,ny/2,24])
    elcX[count] = q[nx/2+1,ny/2,22]
    elcO[count] = q[0,ny/2,22]
    psiCont[count] = dx*sum(psiDiff[0:nx/2])

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

plot(tm, elcX-elcO, 'r-', label='Ez')
plot(tm, -psiCont, 'b-', label='div(B) error')
plot(tm, elcX-elcO-psiCont, 'g-', label='Maxwell')
plot(tmDiff, flxDiff, '-ko', label='d\psi/dt')
legend(loc='upper left')
title('$\Delta Ez$ and $d\psi/dt$')
xlabel('Time')
ylabel('Ez')

show()
