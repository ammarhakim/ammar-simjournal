from pylab import *
import tables
import math
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-p', '--plot', action = 'store',
                  dest = 'fileName',
                  help = 'Hdf5 file to plot')

(options, args) = parser.parse_args()
fileName = options.fileName

fh = tables.openFile(fileName)

q = fh.root.StructGridField
rho = q[:,:,:,0]

grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx = lx/nx
dy = ly/ny
dz = lz/nz

# extend rho array to put in "ghost cells"
rhoExt = zeros((rho.shape[0]+2,rho.shape[1]+2,rho.shape[2]+2), float)
rhoExt[1:-1,1:-1,1:-1] = rho

# apply copy BCs: lower faces
rhoExt[0,:,:] = rhoExt[1,:,:]
rhoExt[:,0,:] = rhoExt[:,1,:]
rhoExt[:,:,0] = rhoExt[:,:,1]
# apply copy BCs: upper faces
rhoExt[-1,:,:] = rhoExt[-2,:,:]
rhoExt[:,-1,:] = rhoExt[:,-2,:]
rhoExt[:,:,-1] = rhoExt[:,:,-2]

# compute gradients
gradRhoX = rhoExt[2:,1:-1,1:-1]-rhoExt[0:-2,1:-1,1:-1]
gradRhoY = rhoExt[1:-1,2:,1:-1]-rhoExt[1:-1,0:-2,1:-1]
gradRhoZ = rhoExt[1:-1,1:-1,2:]-rhoExt[1:-1,1:-1,0:-2]

# data for schlieren image
rhoS = sqrt(gradRhoX**2/(2*dx)**2 + gradRhoY**2/(2*dy)**2 + gradRhoZ**2/(2*dz)**2)

# create a new HDF5 file with this data
fhOut = tables.openFile(fileName[:-3]+"_gradRho.h5", "w")

# create grid group, and copy its attributes
fhOut.create_group(fhOut.root, "StructGrid")
fh.root.StructGrid._v_attrs._f_copy(fhOut.root.StructGrid)

# create time group and copy its attributes
fhOut.create_group(fhOut.root, "timeData")
fh.root.timeData._v_attrs._f_copy(fhOut.root.timeData)

# create new array
fhOut.create_array(fhOut.root, "gradRho", rhoS)
# copy over attributes
fh.root.StructGridField._v_attrs._f_copy(fhOut.root.gradRho)

# close file, saving if
fhOut.close()


