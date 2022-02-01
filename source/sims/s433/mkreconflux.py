from pylab import *
import tables

def getMeshGrid(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
    Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

    return meshgrid(X, Y)

def calcReconFlux(fh, XX, YY, Bx):
    tm = fh.root.timeData._v_attrs.vsTime
    Valf = 0.1
    Lx = 4*pi*5.0
    tmAlf = tm/(Lx/Valf)    
    nx, ny = Bx.shape[0], Bx.shape[1]
    dy = 0.5*Lx/ny
    Byflux = dy*sum(abs(Bx[nx/2,:]))
    return tmAlf, Byflux
    

tmPoints = []
reconFlux = []
for i in range(0,51):
    print ("Working on %d .." % i)
    fh = tables.openFile("../s433/s433-is-coal_q_%d.h5" % i)
    Bx = fh.root.StructGridField[:,:,23]
    X, Y = getMeshGrid(fh.root.StructGrid)
    t, recon = calcReconFlux(fh, X, Y, Bx)
    tmPoints.append(t)
    reconFlux.append(recon)
    fh.close()

plot(tmPoints, reconFlux)
savefig("s433-reconFlux.png")
array(tmPoints).tofile("s433-recon-tmPoints", " ")
array(reconFlux).tofile("s433-recon-reconFlux", " ")
show()
