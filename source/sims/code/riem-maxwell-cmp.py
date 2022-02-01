import numpy
import pylab
import tables
import optparse
import math

pylab.rc('text', usetex=True)

def mkPlot(frame):
    fhW = tables.openFile("../s61/s61-riem-wave_q_%d.h5" % frame)
    fhF_E = tables.openFile("../s62/s62-riem-fdtd_electricField_%d.h5" % frame)
    fhF_B = tables.openFile("../s62/s62-riem-fdtd_magneticField_%d.h5" % frame)

    gridGrp = fhW.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    Xc = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
    Xe = pylab.linspace(xlo, xup-dx, nx)

    q = fhW.root.StructGridField
    F_E = fhF_E.root.StructGridField
    F_B = fhF_B.root.StructGridField

    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.25)
    fig.subplots_adjust(wspace=0.25)

    pylab.subplot(2, 2, 1)
    pylab.plot(Xc, q[:,1], '-k')
    pylab.plot(Xc, F_E[:,1], '-m')
    pylab.title(r"$E_y$")

    pylab.subplot(2, 2, 2)
    pylab.plot(Xc, q[:,2], '-k')
    pylab.plot(Xe, F_E[:,2], '-m')
    pylab.title(r"$E_z$")    

    pylab.subplot(2, 2, 3)
    pylab.plot(Xc, q[:,4], '-k')
    pylab.plot(Xe, F_B[:,1], '-m') 
    pylab.title(r"$B_y$")

    pylab.subplot(2, 2, 4)
    pylab.plot(Xc, q[:,5], '-k')
    pylab.plot(Xc, F_B[:,2], '-m')
    pylab.title(r"$B_z$")

    pylab.savefig("riem-maxwell-cmp_%d.png" % frame)
    
def main():
    mkPlot(2)
    pylab.show()

if __name__ == '__main__': main()
