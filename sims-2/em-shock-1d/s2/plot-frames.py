from pylab import *
import gkedata
import gkedgbasis

for i in range(0,51):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s2-em-shock_distfElc_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm3DPolyOrder1Basis(d)
    X, Y, Z, fv = dg.project(0)
    nx, nvx, nvy = X.shape[0], X.shape[1], X.shape[2]
    figure(1)
    subplot(2,1,1)
    pcolormesh(transpose(fv[:,:,nvy/2]))
    axis('tight')
    colorbar()
    subplot(2,1,2)
    pcolormesh(transpose(fv[:,nvx/2,:]))
    axis('tight')
    colorbar()    
    savefig('s2-em-shock_distfElc_X_VX_%05d.png' % i)
    close()
    d.close()

for i in range(0,51):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s2-em-shock_distfIon_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm3DPolyOrder1Basis(d)
    X, Y, Z, fv = dg.project(0)
    nx, nvx, nvy = X.shape[0], X.shape[1], X.shape[2]
    figure(1)
    subplot(2,1,1)
    pcolormesh(transpose(fv[:,:,nvy/2]))
    axis('tight')
    colorbar()
    subplot(2,1,2)
    pcolormesh(transpose(fv[:,nvx/2,:]))
    axis('tight')
    colorbar()    
    savefig('s2-em-shock_distfIon_X_VX_%05d.png' % i)
    close()
    d.close()    
