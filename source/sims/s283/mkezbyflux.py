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

start = 0
end  = 24
nFrame = end-start+1
tm = zeros((nFrame,), float)
elcX = zeros((nFrame,), float)
elcO = zeros((nFrame,), float)
flx = zeros((nFrame,), float)
psiCont = zeros((nFrame,), float)

count = 0
for i in range(start, end+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("s283-gem-tenmom_q_%d.h5" % i)
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]
    YI = 10

    X = linspace(-LX/2, LX/2, nx)
    Y = linspace(-LY/2, LY/2, ny)

    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    
    psiB_Up = q[:,YI+1,27]
    psiB_Dn = q[:,YI-1,27]
    psiDiff = (psiB_Up-psiB_Dn)/(2*dy)

    tm[count] = fh.root.timeData._v_attrs['vsTime']
    flx[count] = dx*sum(q[1:nx/2,YI,24])
    elcX[count] = 0.5*(q[nx/2,YI,22]+q[nx/2+1,YI,22])
    elcO[count] = 0.5*(q[0,YI,22]+q[1,YI,22])
    psiCont[count] = dx*sum(psiDiff[1:nx/2])

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

plot(tm, elcX, 'r-', label='EzX')
plot(tm, -elcO, 'y-', label='EzO')
plot(tm, -psiCont, 'b-', label='-div(B) error')
plot(tm, elcX-elcO-psiCont, 'g-', label='Total')
plot(tmDiff, flxDiff, '-ko', label='d\psi/dt')
legend(loc='upper right')
title('$\Delta Ez$ and $d\psi/dt$')
xlabel('Time')
ylabel('Ez')

show()
