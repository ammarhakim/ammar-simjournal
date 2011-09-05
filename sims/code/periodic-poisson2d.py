import numpy
import pylab
import tables
import optparse

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Base name of Lua program")

    options, args = parser.parse_args()
    srcFileNm = options.input + "_src.h5"
    solFileNm = options.input + "_sol.h5"

    # open HDF5 files for reading
    srcFile = tables.openFile(srcFileNm)
    solFile = tables.openFile(solFileNm)

    # read grid information
    gridGrp = srcFile.root.StructGrid
    xlo, ylo = gridGrp._v_attrs.vsLowerBounds
    xup, yup = gridGrp._v_attrs.vsUpperBounds
    nx, ny = gridGrp._v_attrs.vsNumCells
    dx = (xup-xlo)/nx
    dy = (yup-ylo)/ny
    X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)
    Y = pylab.linspace(ylo+0.5*dy, yup-0.5*dy, ny)
    XX, YY = pylab.meshgrid(X, Y)

    # read data
    src = srcFile.root.StructGridField
    sol = solFile.root.StructGridField    
    pylab.figure(1)
    pylab.subplot(1,2,1)
    pylab.pcolormesh(XX, YY, src[:,:,0]); pylab.axis('image')
    pylab.subplot(1,2,2)
    pylab.pcolormesh(XX, YY, sol[:,:,0]); pylab.axis('image')
    pylab.contour(XX, YY, sol[:,:,0], 20, colors='k', linestyles='solid')
    pylab.savefig(options.input + "_2d_src_sol.png")

    pylab.show()

if __name__ == '__main__': main()
