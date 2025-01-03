import numpy
import pylab
import tables
import optparse

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Base name of Lua program")
    parser.add_option("-e", "--exact", dest="exact", help="Base name of exact solution input file")
    parser.add_option("-g", "--gamma", dest="gamma", help="Gas adiabatic constant")

    options, args = parser.parse_args()
    consFileNm = options.input + "_q_1.h5"

    gasGamma = float(options.gamma)

    # open HDF5 file for reading
    consFile = tables.openFile(consFileNm)

    # read grid information
    gridGrp = consFile.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    # read data
    cons = consFile.root.StructGridField.read()
    ex_density = pylab.loadtxt(options.exact + "-density.txt")
    ex_velocity = pylab.loadtxt(options.exact + "-velocity.txt")
    ex_pressure = pylab.loadtxt(options.exact + "-pressure.txt")
    ex_ie = pylab.loadtxt(options.exact + "-internal-energy.txt")

    pylab.figure(1)
    pylab.subplot(2, 2, 1)
    pylab.plot(X, cons[:,0], 'k-')
    pylab.plot(ex_density[:,0], ex_density[:,1], 'r-')
    pylab.ylabel("Density")

    pylab.subplot(2, 2, 2)
    pylab.plot(X, cons[:,1]/cons[:,0], 'k-')
    pylab.plot(ex_velocity[:,0], ex_velocity[:,1], 'r-')
    pylab.ylabel("Velocity")

    pr = (gasGamma-1)*(cons[:,2] - 0.5*cons[:,1]*cons[:,1]/cons[:,0])
    pylab.subplot(2, 2, 3)
    pylab.plot(X, pr, 'k-')
    pylab.plot(ex_pressure[:,0], ex_pressure[:,1], 'r-')
    pylab.ylabel("Pressure")

    pylab.subplot(2, 2, 4)
    pylab.plot(X, pr/cons[:,0]/(gasGamma-1), 'k-')
    pylab.plot(ex_ie[:,0], ex_ie[:,1], 'r-')
    pylab.ylabel("Internal Energy")

    pylab.savefig(options.input + "_exact_cmp.png")

    pylab.show()

if __name__ == '__main__': main()    

    
