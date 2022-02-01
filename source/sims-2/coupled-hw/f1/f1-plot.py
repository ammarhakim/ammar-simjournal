from pylab import *

import postgkyl as pg
data = pg.data.GData('f1-hw-lr_chi_100.h5')
interp = pg.data.GInterpModal(data, 2, 'ns')
iGrid, iValues = interp.interpolate()

pcolormesh(iGrid[0], iGrid[1], transpose(iValues[:,:,0]))
axis('image')
show()


