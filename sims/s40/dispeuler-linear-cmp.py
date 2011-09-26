import numpy
import pylab
import tables
import optparse
import math

gasGamma = 2.0
rho0 = 1.0
p0 = 1.0
U0 = 1e-8
nModes = 9
cs0 = math.sqrt(gasGamma*p0/rho0)
Lambda = 10.0

def exactSol(X, t, Lambda):
    Bz = 1.0
    wc = Lambda*Bz
    u1 = 0*X
    U0 = 1e-8
    for n in range(10):
        kn = 2*math.pi*(2*n+1)
        wn = math.sqrt(kn*kn*cs0*cs0 + wc*wc)

        u1 = u1 - U0/(2*n+1)*numpy.sin(kn*X + wn*t)

    return u1

def makePlot(fh, pnx, pny, pn):
    # create grid
    gridGrp = fh.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    # make plot
    q = fh.root.StructGridField
    pylab.subplot(pnx, pny, pn)
    pylab.plot(X, q[:,1]/q[:,0], 'r-')

    dxExact = (xup-xlo)/1000
    Xex = pylab.linspace(xlo+0.5*dxExact, xup-0.5*dxExact, 1000)
    uEx = exactSol(Xex, 3.0, Lambda)
    # plot exact solution
    pylab.plot(Xex, uEx, '-k')    

def main():

    # open HDF5 files for reading
    fh100 = tables.openFile("s40-dispersive-euler_q_1.h5")
    fh200 = tables.openFile("../s42/s42-dispersive-euler_q_1.h5")
    fh300 = tables.openFile("../s43/s43-dispersive-euler_q_1.h5")
    fh400 = tables.openFile("../s44/s44-dispersive-euler_q_1.h5")

    fig = pylab.figure(1)
    makePlot(fh100, 2, 2, 1)
    makePlot(fh200, 2, 2, 2)
    makePlot(fh300, 2, 2, 3)
    makePlot(fh400, 2, 2, 4)

    pylab.savefig("s40424344-dispeuler-cmp.png")

    pylab.show()

if __name__ == '__main__': main()

