from pylab import *
import tables
import math

nx = ny = 768
X = linspace(-12.5, 12.5, nx)
Y = linspace(-12.5, 12.5, ny)

Ez = []
for i in range(61):
    fh = tables.openFile("s238-gemguide-5m_q_%d.h5" % i)
    Ez = fh.root.StructGridField[:,:,12]
    figure(i)
    plot(X, Ez[768/2,:], '-k')
    savefig('s238-Ez-y_%05d.png' % i)
    close()
