import numpy
from pylab import *
import tables

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

fh = tables.openFile("s279-gem-tenmom_q_30.h5")
tm = fh.root.timeData._v_attrs['vsTime']
q = fh.root.StructGridField
nx, ny = q.shape[0], q.shape[1]

X = linspace(-LX/2, LX/2, nx)
Y = linspace(-LY/2, LY/2, ny)
XX, YY = meshgrid(X, Y)

ez = q[:,:,22]
ne = q[:,:,0]/elcMass
ni = q[:,:,10]/ionMass

rhoc = elcCharge*ne+ionCharge*ni
def calcDiv(X, Y, fld):
    dx = X[1]-X[0]
    dy = Y[1]-Y[0]
    divFld = 0.0*fld
    for i in range(1,fld.shape[0]-1):
        for j in range(1,fld.shape[1]-1):
            divFld[i,j] = (fld[i+1,j]-fld[i-1,j])/dx + (fld[i,j+1]-fld[i,j-1])/dy
    return divFld

figure(1)
divE = calcDiv(X, Y, ez)
pcolormesh(XX, YY, transpose(divE-rhoc))
colorbar()

show()
