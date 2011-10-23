import pylab
import tables
import math
import numpy
#pylab.rc('text', usetex=True)

fh = tables.openFile("../s70/s70-plasmabeach_q_100.h5")
q800 = fh.root.StructGridField
dx = 1/800.
X = pylab.linspace(0+0.5*dx, 1-0.5*dx, 800)
pylab.plot(X, q800[:,11], 'm-')

fh = tables.openFile("s69-plasmabeach_q_100.h5")
q400 = fh.root.StructGridField
dx = 1/400.
X = pylab.linspace(0+0.5*dx, 1-0.5*dx, 400)
pylab.plot(X, q400[:,11], 'r-')

fh = tables.openFile("../s68/s68-plasmabeach_q_100.h5")
q200 = fh.root.StructGridField
dx = 1/200.
X = pylab.linspace(0+0.5*dx, 1-0.5*dx, 200)
pylab.plot(X, q200[:,11], 'k-')

fh = tables.openFile("../s67/s67-plasmabeach_q_100.h5")
q100 = fh.root.StructGridField
dx = 1/100.
X = pylab.linspace(0+0.5*dx, 1-0.5*dx, 100)
pylab.plot(X, q100[:,11], 'b-')

pylab.xlabel('X')
pylab.ylabel(r'Ey')

pylab.savefig("plasmabeach_Ey_cmp.png")

pylab.show()


