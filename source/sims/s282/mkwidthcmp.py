import numpy
from pylab import *
import tables

rc('text', usetex=True)

def findXloc(X, fx, val):
    cross = []
    for i in range(X.shape[0]-1):
        if val>fx[i] and val<fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))
        elif val<fx[i] and val>fx[i+1]:
            cross.append(0.5*(X[i]+X[i+1]))            
    return cross

LX = 50.0
LY = 25.0
B0 = 1/15.0
n0 = 1.0
elcCharge = -1.0
ionCharge = 1.0
mu0 = 1.0
ionMass = 1.0
elcMass = ionMass/25
va = B0/sqrt(mu0*elcMass*n0)
ionCycl = ionCharge*B0/ionMass

fh = tables.openFile("s282-gem-tenmom_q_30.h5")
tm = fh.root.timeData._v_attrs['vsTime']*ionCycl
q10 = fh.root.StructGridField
nx, ny = q10.shape[0], q10.shape[1]

fh = tables.openFile("../s280/s280-gemguide-5m_q_30.h5")
q5 = fh.root.StructGridField

X = linspace(-LX/2, LX/2, nx)
Y = linspace(-LY/2, LY/2, ny)
XX, YY = meshgrid(X, Y)

jz10 = q10[:,:,3]
jz5 = q5[:,:,3]
jz10HalfMax = jz10[jz10[:,ny/2].argmax(), ny/2]/2.0
jz5HalfMax = jz5[jz5[:,ny/2].argmax(), ny/2]/2.0

plot(Y, jz10[nx/2,:], '-r', label='10M')
cross = findXloc(Y, jz10[nx/2,:], jz10HalfMax)
plot([cross[0],cross[1]], [jz10HalfMax, jz10HalfMax], '--k')
text(X[nx/2]+0.5*(cross[1]-cross[0]), jz10HalfMax, '%g' % (cross[1]-cross[0]))

plot(Y, jz5[nx/2,:], '-k', label='5M')
cross = findXloc(Y, jz5[nx/2,:], jz5HalfMax)
plot([cross[0],cross[1]], [jz5HalfMax, jz5HalfMax], '--k')
text(X[nx/2]+0.5*(cross[1]-cross[0]), jz5HalfMax, '%g' % (cross[1]-cross[0]))

axis('tight')
gca().set_xlim([-2,2])

legend()
title ('Electron current profile (transverse) at $t\Omega_{ci}$=%g' % tm )
savefig('s282-s280-5m10m_width_cmp.png')
show()



