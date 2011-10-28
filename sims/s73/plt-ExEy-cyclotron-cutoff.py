import pylab
import tables
import math
import numpy
pylab.rc('text', usetex=True)

xlo = 0.0
xup = 0.14
nx = 400
dx = (xup-xlo)/nx
X = pylab.linspace(xlo+0.5*dx, xup-0.5*dx, nx)

driveF = 15.0e9
tEnd = 50.0e-9
nDump = 4
T = pylab.linspace(0, tEnd, nDump+1)


for i in [4]:
    fh = tables.openFile("s73-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField
    
    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.25)

    pylab.subplot(2, 1, 1)
    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    ylims = currAx.get_ylim()
    #currAx.set_ylim([-0.05, 0.05])
    pylab.plot([0.04, 0.04], [ylims[0], ylims[1]], 'k--')
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_x$')

    pylab.subplot(2, 1, 2)
    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-0.05, 0.05])
    pylab.plot([0.04, 0.04], [-0.05, 0.05], 'k--')
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_x$')    

    pylab.savefig('s73-Ex.png')

for i in [4]:
    fh = tables.openFile("s73-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField

    fig = pylab.figure(2)
    
    pylab.plot(X, q[:,11], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-1, 1])
    pylab.plot([0.04, 0.04], [-1, 1], 'k--')
    pylab.axis('tight')
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_y$')    

    pylab.savefig('s73-Ey.png')

for i in [4]:
    fh = tables.openFile("s73-cyclotron-cutoff_q_%d.h5" % i)
    q = fh.root.StructGridField
    
    fig = pylab.figure(3)

    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    ylims = currAx.get_ylim()
    #currAx.set_ylim([-0.05, 0.05])
    pylab.plot([0.04, 0.04], [ylims[0], ylims[1]], 'k--')
    pylab.title("t=%g" % T[i])
    pylab.ylabel(r'$E_x$')

    inset = pylab.axes([0.45, 0.5, 0.4, 0.3])
    pylab.plot(X, q[:,10], 'm-')
    currAx = pylab.gca()
    currAx.set_ylim([-0.05, 0.05])
    pylab.plot([0.04, 0.04], [-0.05, 0.05], 'k--')
    #pylab.ylabel(r'$E_x$')    

    pylab.savefig('s73-Ex-inset.png')
    
pylab.show()
pylab.close()
