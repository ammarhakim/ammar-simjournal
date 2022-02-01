import numpy
import pylab
import tables
import optparse
import math

def mkPlot(frame):

    fhF = tables.openFile("../s64/s64-pulsebox-fdtd-dual_electricField_%d.h5" % frame)
    fhE = tables.openFile("../s59/s59-pulsebox-wave_q_%d.h5" % frame)

    # FDTD grid
    gridGrp = fhF.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    ny = gridGrp._v_attrs.vsNumCells[1]

    q = fhF.root.StructGridField
    # plot FDTD solution
    EzMid = q[:,ny/2,2]
    pylab.plot(X, EzMid, 'ko')

    # wave-grid
    gridGrp = fhE.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    ny = gridGrp._v_attrs.vsNumCells[1]

    # plot wave-solution
    q = fhE.root.StructGridField
    EzMid = 0.5*(q[:,ny/2-1,2] + q[:,ny/2,2])
    pylab.plot(X, EzMid, 'r-')

def main():

    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.25)
    fig.subplots_adjust(wspace=0.25)

    pylab.subplot(2, 1, 1)
    mkPlot(1)
    pylab.title("Electric field at t=1.5")
    pylab.subplot(2, 1, 2)
    mkPlot(2)
    pylab.title("Electric field at t=3.0")    

    pylab.savefig("pulsebox-dual-cmp_1.png")
    
    pylab.show()

if __name__ == '__main__': main()

