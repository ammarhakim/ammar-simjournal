import tables
import pylab
import math
import numpy

baseName = "s212-advect-diffuse"
alpha = 0.1
tend = 2*math.pi

fh = tables.openFile(baseName+"_q_0.h5")
q0 = fh.root.StructGridField

# read in grid data
grid = fh.root.StructGrid
lower = grid._v_attrs.vsLowerBounds
upper = grid._v_attrs.vsUpperBounds
cells = grid._v_attrs.vsNumCells

nx = cells[0]
dx = (upper[0]-lower[0])/cells[0]
Xlo = pylab.linspace(lower[0], upper[0]-dx, nx)
X = pylab.linspace(lower[0], upper[0], nx+1)

fh = tables.openFile(baseName+"_q_4.h5")
q4 = fh.root.StructGridField

def exactSol(X, tend):
    return pylab.exp(-alpha*tend)*pylab.sin(X-tend)

qEx = exactSol(Xlo, tend)

error = numpy.abs(qEx-q4[:,0,0]).sum()/nx
print dx, error

def drawLines(X, q, color):
    nx = X.shape[0]
    for i in range(nx-1):
        pylab.plot([X[i], X[i+1]], [q[i,0,0], q[i,0,1]], color)

pylab.figure(1)
# make plots
pylab.plot(X, exactSol(X, tend), 'ro')
drawLines(X, q4, '-k')
drawLines(X, q0, '-r')
pylab.axis('tight')
pylab.xlabel('X')

pylab.savefig(baseName+"-cmp.png")
pylab.show()
