import numpy
import pylab
import tables
import optparse
import math

def mkPlot(frame):

    fhW = tables.openFile("../s57/s57-pulsebox-wave_q_%d.h5" % frame)
    fhF = tables.openFile("../s58/s58-pulsebox-fdtd_electricField_%d.h5" % frame)

    fhE = tables.openFile("../s59/s59-pulsebox-wave_q_%d.h5" % frame)

    # wave-grid
    gridGrp = fhW.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    ny = gridGrp._v_attrs.vsNumCells[1]

    # plot wave-solution
    q = fhW.root.StructGridField
    EzMid = 0.5*(q[:,ny/2-1,2] + q[:,ny/2,2])
    pylab.plot(X, EzMid, 'k-')

    # FDTD grid
    gridGrp = fhF.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo, xup-dx, nx)

    ny = gridGrp._v_attrs.vsNumCells[1]    

    q = fhF.root.StructGridField
    # plot FDTD solution
    EzMid = q[:,ny/2,2]
    pylab.plot(X, EzMid, 'm-')

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

def mkPlot2d(frame):

    fhW = tables.openFile("../s57/s57-pulsebox-wave_q_%d.h5" % frame)

    # wave-grid
    gridGrp = fhW.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    ylo = gridGrp._v_attrs.vsLowerBounds[1]
    yup = gridGrp._v_attrs.vsUpperBounds[1]
    ny = gridGrp._v_attrs.vsNumCells[1]
    dy = (yup-ylo)/ny
    Y = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)

    XX, YY = pylab.meshgrid(X, Y)
    q = fhW.root.StructGridField
    Ez = q[:,:,2]

    pylab.pcolormesh(XX, YY, Ez.transpose())
    pylab.axis('image')

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

    pylab.savefig("pulse-box-cmp_1.png")

    fig = pylab.figure(2)
    pylab.subplot(1, 2, 1)
    mkPlot2d(1)
    pylab.title("Electric field at t=1.5")
    pylab.subplot(1, 2, 2)
    mkPlot2d(2)
    pylab.title("Electric field at t=3.0")

    pylab.savefig("pulse-box-cmp_2d.png")    
    
    pylab.show()

if __name__ == '__main__': main()

