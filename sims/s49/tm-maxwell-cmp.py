import numpy
import pylab
import tables
import optparse
import math

class Counter:
    def __init__(self):
        self.count = 0

    def getCount(self):
        self.count = self.count + 1
        return self.count

counter = Counter()

def getExact(XX, YY, tm):
    m, n = 8, 5
    Xu, Yu = 80, 40
    a = m*math.pi/Xu
    b = n*math.pi/Yu
    c0 = 2.99792458e8
    w = c0*math.sqrt(a*a + b*b)

    return numpy.sin(a*XX)*numpy.sin(b*YY)*math.cos(w*tm)

def mkPlot(luaFl, frame, tm):
    fh = tables.openFile("%s_q_%d.h5" % (luaFl, frame))

    gridGrp = fh.root.StructGrid
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

    q = fh.root.StructGridField
    # make 2D plot of solution
    pylab.figure(counter.getCount())
    pylab.pcolormesh(XX, YY, q[:,:,2].transpose())
    pylab.axis('image')
    pylab.xlabel("X")
    pylab.ylabel("Y")
    pylab.title("Time t=%g" % tm)
    pylab.colorbar()
    pylab.savefig("%s_2d_%d.png" % (luaFl, frame))

    Ez = getExact(XX, YY, tm)
    pylab.figure(counter.getCount())
    pylab.plot(X, q[:,ny/2,2], 'k-')
    pylab.plot(X, Ez.transpose()[:,ny/2], 'r-')
    pylab.axis('tight')
    pylab.savefig("%s_1d_%d.png" % (luaFl, frame))

def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Base name of Lua program")

    # open files for plotting
    options, args = parser.parse_args()
    luaFl = options.input

    mkPlot(luaFl, 1, 75e-9)
    mkPlot(luaFl, 2, 150e-9)

if __name__ == '__main__': main()
