from pylab import *
import gkedata
import gkedgbasis

for i in range(0,51):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s1-em-shock_distfElc_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm2DPolyOrder2Basis(d)
    X, Y, fv = dg.project(0)
    figure(1)
    pcolormesh(transpose(fv))
    axis('tight')
    colorbar()
    savefig('s1-em-shock_distfElc_X_VX_%05d.png' % i)
    close()
    d.close()

for i in range(0,51):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s1-em-shock_distfIon_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm3DPolyOrder2Basis(d)
    X, Y, fv = dg.project(0)
    figure(1)
    pcolormesh(transpose(fv[:,:,nvy/2]))
    axis('tight')
    colorbar()
    savefig('s1-em-shock_distfIon_X_VX_%05d.png' % i)
    close()
    d.close()    
