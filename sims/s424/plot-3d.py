from pylab import *
import tables
import numpy

# customization for figure
#rcParams['lines.linewidth']            = 2
rcParams['font.size']                  = 24
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

def getVars(q):
    gasGamma = 1.4
    rho = q[:,:,:,0]
    u = q[:,:,:,1]/rho
    v = q[:,:,:,2]/rho
    w = q[:,:,:,3]/rho
    Er = q[:,:,:,4]
    ke = 0.5*rho*(u*u+v*v+w*w)
    pr = (gasGamma-1)*(Er - ke)
    cs = sqrt(gasGamma*pr/rho)
    mach = sqrt(u*u+v*v+w*w)/cs
    return rho, pr, mach

# 3D
fh = tables.openFile("s424-super-euler-3d_inOut.h5")
inOut = fh.root.StructGridField
grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx = lx/nx
dy = ly/ny
dz = lz/nz
Xe = linspace(0.0, lx, nx+1)
Ze = linspace(0.0, lz, nz+1)

Xc = linspace(0.5*dx, lx-0.5*dx, nx)
Zc = linspace(0.5*dz, lz-0.5*dz, nz)

fh = tables.openFile("s424-super-euler-3d_q_10.h5")
q = fh.root.StructGridField
rho, pr, mach = getVars(q)

inOutA = inOut[:,0,:]
rhoMa = numpy.ma.masked_where(inOut[:,:,0]<0, rho[:,0,:])

pcolormesh(Xc, Zc, transpose(rhoMa))

axis('image')

savefig('s424-rho-3d.png', dpi=300, bbox_inches='tight')
show()
