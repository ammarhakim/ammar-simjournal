
import numpy as np
import matplotlib.pyplot as plt

'''
 Utility functions
'''	
def calc_phi2d(fx, fy, dx=1, dy=1):
	'''
	Calculate phi by integrating dphi = fx*dx + fy*dy, phi[0,0] = 0 (<---> [fx, fy] = grad(phi))
	'''
	ny,nx=fx.shape
	phi=np.zeros((ny,nx))
	for jx in range(1,nx):
		phi[0,jx]=phi[0,jx-1]+fx[0,jx]*dx
	for jy in range(1,ny):
		phi[jy,:]=phi[jy-1,:]+fy[jy,:]*dy
	return phi

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
		#psi[jy,:]=psi[jy-1,:]+0.5*(fx[jy-1,:]+fx[jy,:])*dy
		psi[jy,:]=psi[jy-1,:]+fx[jy,:]*dy
	# since f = rot(A) gives extra restraints on f (e.g., div(f)=0)
	# it makes sense that information provided by fy[1:nx,:] is useless here
	return psi

def calc_psi2d_Kai(fx,fy, dx=1, dy=1):
	ny,nx=fx.shape
	psi=np.zeros((ny,nx))
	for i in range(1,ny):
	    psi[i,0] = psi[i-1,0] + dy * ( fx[i,0] + fx[i-1,0] )/2.
	for j in range(1,nx):
	    psi[:,j] = psi[:,j-1] - dx * ( fy[:,j-1] + fy[:,j] )/2.
	return psi	
