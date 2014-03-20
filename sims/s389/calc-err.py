from pylab import *
import tables

fh = tables.openFile("s389-euler-ds-2d_q_0.h5")
q0 = fh.root.StructGridField
fh = tables.openFile("s389-euler-ds-2d_q_1.h5")
q1 = fh.root.StructGridField
nx = fh.root.StructGrid._v_attrs.vsNumCells[0]
ny = fh.root.StructGrid._v_attrs.vsNumCells[1]
dx = 2.0/nx
err = sum( (q0[:,:,0]-q1[:,:,0])**2 )/(nx*ny)

print dx, err
