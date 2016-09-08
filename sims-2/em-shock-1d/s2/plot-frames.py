from pylab import *
import gkedata
import gkedgbasis

for i in range(0,41):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s2-em-shock_distfElc_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm3DPolyOrder2Basis(d)
    X, Y, Z, fv = dg.project(0)
    pcolormesh(transpose(fv[:,24,:]))
    axis('tight')
    colorbar()
    savefig('s2-em-shock_distfElc_X_VX_%05d.png' % i)
    close()
    d.close()

for i in range(0,41):
    print("Working on %d ..." % i)
    d = gkedata.GkeData("s2-em-shock_distfIon_%d.h5" % i)
    dg = gkedgbasis.GkeDgSerendipNorm3DPolyOrder2Basis(d)
    X, Y, Z, fv = dg.project(0)
    pcolormesh(transpose(fv[:,24,:]))
    axis('tight')
    colorbar()
    savefig('s2-em-shock_distfIon_X_VX_%05d.png' % i)
    close()
    d.close()    
