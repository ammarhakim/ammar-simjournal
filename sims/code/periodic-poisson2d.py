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
    solCDFileNm = options.input + "_solCD.h5"

    # open HDF5 files for reading
    srcFile = tables.openFile(srcFileNm)
    solFile = tables.openFile(solFileNm)
    solCDFile = tables.openFile(solCDFileNm)

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

    # for plotting we need vertex coordinates
    Xp = pylab.linspace(xlo, xup, nx+1)
    Yp = pylab.linspace(ylo, yup, ny+1)
    XXp, YYp = pylab.meshgrid(Xp, Yp)

    # read data
    src = srcFile.root.StructGridField
    sol = solFile.root.StructGridField
    solCD = solCDFile.root.StructGridField

    pylab.figure(1)
    pylab.subplot(1,2,1)
    pylab.pcolormesh(XXp, YYp, src[:,:,0].transpose());
    pylab.axis('image')
    pylab.subplot(1,2,2)
    pylab.pcolormesh(XXp, YYp, sol[:,:,0].transpose());
    pylab.axis('image')
    pylab.contour(XX, YY, sol[:,:,0].transpose(), 20, colors='k', linestyles='solid')
    pylab.savefig(options.input + "_2d_src_sol.png")

    src = src.read() - src.read().sum()/(nx*ny) # read adjust source
    pylab.figure(2)
    pylab.subplot(1,2,1)
    pylab.plot(X, solCD[:,ny/2], '-k', X, src[:,ny/2], 'ro')
    pylab.xlabel('X')
    pylab.subplot(1,2,2)
    pylab.plot(Y, solCD[nx/2,:], '-k', Y, src[nx/2,:], 'ro')
    pylab.xlabel('Y')    
    pylab.savefig(options.input + "_1d_CD_cmp.png")

    pylab.show()

if __name__ == '__main__': main()
