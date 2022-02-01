import tables
from pylab import *

phi = []
for i in range(101):
    fh = tables.open_file("s1-es-iaw_em_%d.h5" % i)
    q = fh.root.StructGridField
    phi.append(q[10,0])

plot(array(phi))
show()
