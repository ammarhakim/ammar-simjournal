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

def getExact(X, Y, tm):
    m, n = 8, 5
    Xu, Yu = 80, 40
    a = m*math.pi/Xu
    b = n*math.pi/Yu
    c0 = 2.99792458e8
    w = c0*math.sqrt(a*a + b*b)

    Ez = numpy.zeros( (X.shape[0], Y.shape[0]), numpy.float )
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            Ez[i,j] = math.sin(a*X[i])*math.sin(b*Y[j])*math.cos(w*tm)

    return Ez

def getExactX(X, Y, tm):
    m, n = 8, 5
    Xu, Yu = 80, 40
    a = m*math.pi/Xu
    b = n*math.pi/Yu
    c0 = 2.99792458e8
    w = c0*math.sqrt(a*a + b*b)

    Ez = numpy.sin(a*X)*math.sin(b*Y)*math.cos(w*tm)

    return Ez

flNms = [ ("../s49/s49-tm-maxwell-wave", "../s53/s53-tm-maxwell-fdtd"),
          ("../s50/s50-tm-maxwell-wave", "../s54/s54-tm-maxwell-fdtd"),
          ("../s51/s51-tm-maxwell-wave", "../s55/s55-tm-maxwell-fdtd"),
          ("../s52/s52-tm-maxwell-wave", "../s56/s56-tm-maxwell-fdtd")]

def mkPlot(fPair, frame, tm):

    fhW = tables.openFile("%s_q_%d.h5" % (fPair[0], frame))
    fhF = tables.openFile("%s_electricField_%d.h5" % (fPair[1], frame))

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

    # plot exact solution
    Ez = getExactX(X, 20, tm)    
    pylab.plot(X, Ez, 'r-')

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

def main():

    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.25)
    fig.subplots_adjust(wspace=0.25)

    pylab.subplot(2, 2, 1)
    mkPlot(flNms[0], 1, 75e-9)
    pylab.subplot(2, 2, 2)
    mkPlot(flNms[1], 1, 75e-9)    
    pylab.subplot(2, 2, 3)
    mkPlot(flNms[2], 1, 75e-9)
    pylab.subplot(2, 2, 4)
    mkPlot(flNms[3], 1, 75e-9)

    pylab.suptitle("Electric field at t=75 ns")
    pylab.savefig("tm-maxwell-cmp-1.png")

    fig = pylab.figure(2)
    fig.subplots_adjust(hspace=0.25)
    fig.subplots_adjust(wspace=0.25)

    pylab.subplot(2, 2, 1)
    mkPlot(flNms[0], 2, 150e-9)
    pylab.subplot(2, 2, 2)
    mkPlot(flNms[1], 2, 150e-9)    
    pylab.subplot(2, 2, 3)
    mkPlot(flNms[2], 2, 150e-9)
    pylab.subplot(2, 2, 4)
    mkPlot(flNms[3], 2, 150e-9)

    pylab.suptitle("Electric field at t=150 ns")
    pylab.savefig("tm-maxwell-cmp-2.png")
    
    pylab.show()

if __name__ == '__main__': main()
