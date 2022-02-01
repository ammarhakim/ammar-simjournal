import tables
import pylab
import math
import numpy
from matplotlib import rcParams
import matplotlib.pyplot as plt

# customization for figure
rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 18
#rcParams['xtick.major.size']           = 8 # default is 4
#rcParams['xtick.major.width']          = 3 # default is 0.5
#rcParams['ytick.major.size']           = 8 # default is 4
#rcParams['ytick.major.width']          = 3 # default is 0.5
rcParams['figure.facecolor']           = 'white'
#rcParams['figure.subplot.bottom']      = 0.125
#rcParams['figure.subplot.right']       = 0.85 # keep labels/ticks of colobar in figure
rcParams['image.interpolation']        = 'none'
rcParams['image.origin']               = 'lower'
rcParams['contour.negative_linestyle'] = 'solid'
#rcParams['savefig.bbox']               = 'tight'

def plotLines(X, q0, q1, color):
    for i in range(X.shape[0]-1):
        ql = q0[i] - q1[i]
        qr = q0[i] + q1[i]
        pylab.plot([X[i], X[i+1]], [ql, qr], color)

def calcDeriv(T, val):
    dVal = numpy.zeros((T.shape[0]-1,), numpy.float)
    tVal = numpy.zeros((T.shape[0]-1,), numpy.float) 
    for i in range(T.shape[0]-1):
        tVal[i] = 0.5*(T[i]+T[i+1])
        vAv = 0.5*(val[i+1]+val[i])
        dVal[i] = (val[i+1]-val[i])/(T[i+1]-T[i])
    return tVal, dVal

def exactSol(X, tend):
    return pylab.exp(-alpha*tend)*numpy.sin(X)

def exactSolSS(X):
    return numpy.sin(X)

baseList = ["../s346/s346-modal-dg-diffuse",
            "../s336/s336-modal-dg-diffuse", "../s340/s340-modal-dg-diffuse"]
count = 0

fig = pylab.figure(1)
for baseName in baseList:
    pylab.subplot(3, 1, count+1)

    alpha = 1.0
    tend = 1.0

    fh = tables.openFile(baseName+"_q_0.h5")
    q0 = fh.root.StructGridField

    grid = fh.root.StructGrid
    lower = grid._v_attrs.vsLowerBounds
    upper = grid._v_attrs.vsUpperBounds
    cells = grid._v_attrs.vsNumCells

    nx = cells[0]
    dx = (upper[0]-lower[0])/cells[0]
    X = pylab.linspace(lower[0], upper[0], nx+1)

    fh = tables.openFile(baseName+"_q_1.h5")
    q = fh.root.StructGridField

    plotLines(X, q[:,0], q[:,1], '-r')

    Xhr = pylab.linspace(lower[0], upper[0], 1000)
    qExact = exactSol(Xhr, 3.0/40.0)
    pylab.plot(Xhr, qExact, '-m', linewidth=1.0)
    #pylab.title('Steady-state solution for %s scheme' % titleStr[count] )
    #pylab.xlabel('X')
    #pylab.ylabel('T')
    pylab.axis('tight')
    count = count + 1

pylab.savefig('s346-s336-s340-dg-diffuse.pdf')
pylab.show()


