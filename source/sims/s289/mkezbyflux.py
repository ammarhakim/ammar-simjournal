import numpy
from pylab import *
import tables

rc('text', usetex=True)

LX = 2.0
LY = 2.0

start = 0
end  = 100
nFrame = end-start+1
tm = zeros((nFrame,), float)
elcX = zeros((nFrame,), float)
elcO = zeros((nFrame,), float)
flx = zeros((nFrame,), float)
psiCont = zeros((nFrame,), float)

count = 0
for i in range(start, end+1):
    print ("Working on %d ..." % i)
    fh = tables.openFile("s289-pulse-box-5m_q_%d.h5" % i)
    q = fh.root.StructGridField
    nx, ny = q.shape[0], q.shape[1]
    YI = ny/2

    X = linspace(-LX/2, LX/2, nx)
    Y = linspace(-LY/2, LY/2, ny)

    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    
    psiB_Up = q[:,YI+1,17]
    psiB_Dn = q[:,YI-1,17]
    psiDiff = (psiB_Up-psiB_Dn)/(2*dy)

    tm[count] = fh.root.timeData._v_attrs['vsTime']
    flx[count] = dx*sum(q[1:nx/2,YI,14])
    elcX[count] = 0.5*(q[nx/2,YI,12]+q[nx/2+1,YI,12])
    elcO[count] = 0.5*(q[0,YI,12]+q[1,YI,12])
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
plot(tmDiff, flxDiff, '-k', label='d\psi/dt')
legend(loc='best')
title('$\Delta Ez$ and $d\psi/dt$')
xlabel('Time')
ylabel('Ez')

savefig('s284-ezbyflux.png')

show()
