from mayavi import mlab
from mayavi.mlab import *
from mayavi.mlab import savefig
import pylab
import tables

def getMesh(grid):
    xlo, ylo, zlo = grid._v_attrs.vsLowerBounds
    xup, yup, zup = grid._v_attrs.vsUpperBounds
    nx, ny, nz = grid._v_attrs.vsNumCells
    x, y, z = pylab.mgrid[xlo:xup:nx*1j, ylo:yup:ny*1j, zlo:zup:nz*1j]
    return x, y, z

for i in range(10,11):
    print ("Working on frame %d ... " % i)
    fh = tables.openFile("s2-em-shock_distfElc_%d.h5" % i)
    grid = fh.root.StructGrid
    nx, ny, nz = grid._v_attrs.vsNumCells
    x, y, z = getMesh(grid)

    rho = fh.root.StructGridField[:,:,:,0]
    rhoSd = pipeline.scalar_field(x,y,z,rho)

    # add various figures
    pipeline.image_plane_widget(rhoSd,
                                plane_orientation='x_axes', 
                                slice_index=nx/2)
    pipeline.image_plane_widget(rhoSd,
                                plane_orientation='x_axes', 
                                slice_index=nx/4)
    pipeline.image_plane_widget(rhoSd,
                                plane_orientation='x_axes', 
                                slice_index=3*nx/4)

    outline()
    colorbar(orientation='vertical')
    pipeline.iso_surface(rhoSd, contours=[6.0], opacity=0.75)
    roll(0.0)

    savefig('s446-euler-rt_rho_%05d.png' % i, magnification=2.0)
    close()
    fh.close()
