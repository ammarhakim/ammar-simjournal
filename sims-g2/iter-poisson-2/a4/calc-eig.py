import scipy.io
import scipy.sparse.linalg

stiffMat = scipy.io.mmread("a3-direct-poisson-32-stiffMat.mm")
vals_X, vecs_X = scipy.sparse.linalg.eigs(stiffMat, which='LR')
vals_N, vecs_N = scipy.sparse.linalg.eigs(stiffMat, which='SR')

lam_max = vals_X[1].real # skip the positive eigenvalue
lam_min = vals_N[5].real

print( "L_max %g. L_min %g" % (lam_max, lam_min) )
