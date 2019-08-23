from pylab import *
import postgkyl as pg
import numpy

style.use('postgkyl.mplstyle')

baseNm = "n1-es-buneman"

def getJ(fr):
    data = pg.GData("%s_elc_M1i_%d.bp" % (baseNm, fr))
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, m1e = dg.interpolate()
    
    data = pg.GData("%s_ion_M1i_%d.bp" % (baseNm, fr))
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, m1i = dg.interpolate()

    Jg = 0.159*ones(Ji.shape, float)

    Xn = XX[0]; dx = Xn[1]-Xn[0]
    Xc = linspace(Xn[0]+0.5*dx, Xn[-1]-0.5*dx, Xn.shape[0]-1)

    return Xc, m1i, -m1e+Jg
