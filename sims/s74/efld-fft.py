import pylab
import tables
import math
import numpy
import findpeaks

nx = 400
nt = 191
Ex_TX = numpy.zeros((nt+1, nx), numpy.float)

for i in range(191):
    fh = tables.openFile("s74-icw_EM_%d.h5" % i)
    q = fh.root.StructGridField
    Ex = q[:,0]
    Ex_TX[i,:] = Ex
