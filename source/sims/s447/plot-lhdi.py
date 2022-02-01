from pylab import *
import tables

fh = tables.openFile("s447-5m-lhdi_q_0.h5")
nx, ny = fh.root.StructGrid._v_attrs.vsNumCells[0], fh.root.StructGrid._v_attrs.vsNumCells[1]
xl, yl = fh.root.StructGrid._v_attrs.vsLowerBounds[0], fh.root.StructGrid._v_attrs.vsLowerBounds[1]
xu, yu = fh.root.StructGrid._v_attrs.vsUpperBounds[0], fh.root.StructGrid._v_attrs.vsUpperBounds[1]

dx = (xu-xl)/nx
dy = (yu-yl)/ny

X = linspace(xl+0.5*dx, xu-0.5*dx, nx)
Y = linspace(yl+0.5*dy, yu-0.5*dy, ny)

for i in range(0,51):
    print ("Working on frame %d ..." % i)
    fh = tables.openFile("s447-5m-lhdi_q_%d.h5" % i)
    q = fh.root.StructGridField
    figure(1)
    pcolormesh(X, Y, transpose(q[:,:,0]))
    axis('image')
    savefig('s447-5m-lhdi-rho_%05d.png' % i)
    close()
    fh.close()

