from pylab import *
import postgkyl as pg
import numpy
import numpy as np

style.use('postgkyl.mplstyle')

def _colorbar(obj, fig, ax, label=""):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    return fig.colorbar(obj, cax=cax, label=label)

def calc_psi2d(fx, fy, dx=1, dy=1):#solenoidal flows
	'''
	Calcualte psi by integrating dpsi = -fy*dx + fx*dy, psi[0,0]=0.
	Notes: 
		1. (fx=dpsi/dy,fy=-dpsi/dx) is called Hamiltonian gradient of psi, and	contours of psi give vector field (fx, fy);
		2. div(f)=0
	'''
	#FIXME: avoid loop for higher performance
	ny,nx=fx.shape
	psi=np.zeros((ny,nx))
	for jx in range(1,nx):
		#psi[0,jx]=psi[0,jx-1]-0.5*(fy[0,jx-1]+fy[0,jx])*dx
		psi[0,jx]=psi[0,jx-1]-fy[0,jx]*dx
	for jy in range(1,ny):
		psi[jy,:]=psi[jy-1,:]+fx[jy,:]*dy
	# since f = rot(A) gives extra restraints on f (e.g., div(f)=0)
	# it makes sense that information provided by fy[1:nx,:] is useless here
	return psi

for i in range(399, 400):
    
    print("Working on %d ... " % i)

    fig = figure(1)

    # elc
    data = pg.GData("ot1-orsag-tang_elc_M1i_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q0 = dg.interpolate(2)
    X, Y = meshgrid(XX[0], XX[1])
    im = pcolormesh(X, Y, transpose(q0[:,:,0]))

    # psi
    data = pg.GData("ot1-orsag-tang_field_%d.bp" % i)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, Bx = dg.interpolate(3)
    XX, By = dg.interpolate(4)

    X = XX[0]; Y = XX[1]
    dx = X[1]-X[0]; dy = Y[1]-Y[0]

    Xc = linspace(X[0]+0.5*dx, X[-1]-0.5*dx, X.shape[0]-1)
    Yc = linspace(Y[0]+0.5*dy, Y[-1]-0.5*dy, Y.shape[0]-1)
    
    psiB = calc_psi2d(transpose(Bx[:,:,0]), transpose(By[:,:,0]), dx, dy)

    contour(Xc, Yc, psiB, 20, colors='w', linestyles='solid', linewidths=1)
    
    title(r"$\Omega t = {:0.2f}$".format(20*data.time/5000.0))
    axis('image')

    savefig("ot1-orsag-tang_elc_M1i_c2_%05d.png" % i, dpi=200)
    close()
    

