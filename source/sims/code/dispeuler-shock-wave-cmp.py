import numpy
import pylab
import tables
import optparse

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Base name of Lua program")
    parser.add_option("-g", "--gamma", dest="gamma", help="Gas adiabatic constant")

    options, args = parser.parse_args()
    primFileNm = options.input + "_prim_1.h5"

    gasGamma = float(options.gamma)

    # open HDF5 file for reading
    primFile = tables.openFile(primFileNm)

    # read grid information
    gridGrp = primFile.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    prim = primFile.root.StructGridField

    pylab.figure(1)
    pylab.subplot(2, 2, 1)
    pylab.plot(X, prim[:,0], 'k-')
    pylab.ylabel("Density")

    pylab.subplot(2, 2, 2)
    pylab.plot(X, prim[:,1], 'k-')
    pylab.ylabel("Velocity")

    pylab.subplot(2, 2, 3)
    pylab.plot(X, prim[:,4], 'k-')
    pylab.ylabel("Pressure")

    pylab.subplot(2, 2, 4)
    pylab.plot(X, prim[:,4]/prim[:,0]/(gasGamma-1), 'k-')
    pylab.ylabel("Internal Energy")

    pylab.savefig(options.input + "_sol.png")

    pylab.show()

if __name__ == '__main__': main()    

    
