import numpy
import pylab
import tables
import optparse

pylab.rc('text', usetex=True)

def mkPlot(fh100, fh200, tm):
    q100 = fh100.root.StructGridField
    q200 = fh200.root.StructGridField

    dx100 = 1/100.0
    X100 = pylab.linspace(0.5*dx100, 1-0.5*dx100, 100)

    dx200 = 1/200.0
    X200 = pylab.linspace(0.5*dx200, 1-0.5*dx200, 200)

    pylab.plot(X100, q100[:,1], 'r-')
    pylab.plot(X200, q200[:,1], 'k-')
    pylab.xlabel('X')
    pylab.ylabel(r'$E_y$')
    pylab.title('t=%s ns' % tm)

def main():
    fh100 = tables.openFile("s65-plasmabeach-maxwell_q_1.h5")
    fh200 = tables.openFile("../s66/s66-plasmabeach-maxwell_q_1.h5")

    fig = pylab.figure(1)
    fig.subplots_adjust(hspace=0.3)
    #fig.subplots_adjust(wspace=0.5)
    
    pylab.subplot(2, 1, 1)
    mkPlot(fh100, fh200, '2.5')

    fh100 = tables.openFile("s65-plasmabeach-maxwell_q_2.h5")
    fh200 = tables.openFile("../s66/s66-plasmabeach-maxwell_q_2.h5")

    pylab.subplot(2, 1, 2)
    mkPlot(fh100, fh200, '5.0')

    pylab.savefig('plasmabeach-maxwell-cmp.png')
    pylab.show()

if __name__ == '__main__': main()
