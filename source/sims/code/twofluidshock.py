import numpy
import pylab
import tables
import optparse

charge = 1.0
gasGamma = 5./3.
elcCharge = -charge
ionCharge = charge
elcMass = 1/1836.2
ionMass = 1.0
lightSpeed = 1.0
epsilon0 = 1.0

def getPressure(q, s):
    Er = q[:,s+4] - 0.5*(q[:,s+1]**2 + q[:,s+2]**2 + q[:,s+3]**2)/q[:,s+0]
    return Er*(gasGamma-1)

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Base name of Lua program")

    options, args = parser.parse_args()
    fileNm = options.input + "_q_1.h5"

    # open HDF5 file for reading
    fh = tables.openFile(fileNm)

    # create grid
    gridGrp = fh.root.StructGrid
    xlo = gridGrp._v_attrs.vsLowerBounds[0]
    xup = gridGrp._v_attrs.vsUpperBounds[0]
    nx = gridGrp._v_attrs.vsNumCells[0]
    dx = (xup-xlo)/nx
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

    # make plots
    q = fh.root.StructGridField
    ne = q[:,0]/elcMass # electron number density
    ni = q[:,5]/ionMass # ion number density

    pylab.figure(1)
    pylab.plot(X, ne, 'r-', X, ni, 'k-')
    pylab.xlabel("X")
    pylab.ylabel("Number Density")

    pylab.savefig(options.input + "_neni.png")

    fig = pylab.figure(2)
    fig.subplots_adjust(hspace=0.5)
    fig.subplots_adjust(wspace=0.5)
    
    ue = q[:,1]/q[:,0]
    ui = q[:,6]/q[:,5]
    pe = getPressure(q, 0)
    pi = getPressure(q, 5)

    pylab.subplot(2, 2, 1)
    pylab.plot(X, ue, '-r')
    pylab.xlabel("X")
    pylab.ylabel("X Velocity")

    pylab.subplot(2, 2, 2)
    pylab.plot(X, ui, '-k')
    pylab.xlabel("X")
    pylab.ylabel("X Velocity")    

    pylab.subplot(2, 2, 3)
    pylab.plot(X, pe, '-r')
    pylab.xlabel("X")
    pylab.ylabel("Pressure")    

    pylab.subplot(2, 2, 4)
    pylab.plot(X, pi, '-k')
    pylab.xlabel("X")
    pylab.ylabel("Pressure")

    pylab.savefig(options.input + "_up.png")

    pylab.show()

if __name__ == '__main__': main()
