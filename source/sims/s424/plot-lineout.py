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

# function to compute lineout along ray with angle theta
def calcLineOut(nv, dx, dy, rho, theta):
    vals = zeros((nv, 1), float)
    cost = cos(theta)
    sint = sin(theta)
    dv = 1.0/nv
    vp = linspace(0.5*dv, 1.0-0.5*dv, nv)
    xp = cost*vp
    yp = sint*vp
    for i in range(nv):
        ix = (int) (xp[i]/dx)
        iy = (int) (yp[i]/dy)
        vals[i] = rho[ix,iy]
    return vals

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

# axi-symmetric

fh = tables.openFile("../s423/s423-euler-blunt-rz_q_10.h5")
grid = fh.root.StructGrid
nx, ny, nz = grid._v_attrs.vsNumCells[0], grid._v_attrs.vsNumCells[1], grid._v_attrs.vsNumCells[2]
lx, ly, lz = grid._v_attrs.vsUpperBounds[0], grid._v_attrs.vsUpperBounds[1], grid._v_attrs.vsUpperBounds[2]
dx_ax = lx/nx
dy_ax = ly/ny
dz_ax = lz/nz

X_ax = linspace(0, 1.0, nx)
rho_ax, pr_ax, mach_ax = getVars(fh.root.StructGridField)

iz_ax = 1.25/dz_ax

# 3D
fh = tables.openFile("s424-super-euler-3d_inOut.h5")
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

fh = tables.openFile("s424-super-euler-3d_q_%d.h5" % 10)
q = fh.root.StructGridField
rho, pr, mach = getVars(q)

nv = 100
Xr = linspace(0, 1.0, nv)
iz = 1.25/dz

figure(1)
tRange = linspace(0, pi/4, 20)
for i in range(tRange.shape[0]):
    v = calcLineOut(nv, dx, dy, pr[:,:,iz], tRange[i])
    plot(Xr, v, '-b')
    plot(X_ax, pr_ax[:,0,iz_ax], '-k', linewidth=2)

xlabel('Radial Location')
ylabel('Pressure')

savefig('s424-lineout-cmp.png', dpi=300, bbox_inches='tight')
show()
