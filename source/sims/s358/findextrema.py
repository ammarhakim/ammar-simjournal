import numpy
import pylab
import math

pylab.rc('text', usetex=True)

def areSignDiff(a, b):
    if a>0 and b>0:
        return False
    elif a<0 and b<0:
        return False
    return True

def findextrema(fld):
    nx = fld.shape[0]
    idxLst = []
    prevSlope = fld[1]-fld[0]
    for i in range(1,nx-1):
        slope = fld[i+1]-fld[i]
        if areSignDiff(prevSlope, slope):
            idxLst.append(i)
        prevSlope = slope
    return idxLst

def makePlot(T, fld):
    extrem = findextrema(numpy.log(fld))
    pylab.semilogy(T, fld, 'b-')
    ylim = pylab.gca().get_ylim()
    #for e in extrem:
    #    pylab.plot([T[e],T[e]], [ylim[0], ylim[1]], 'k-')
    exArr = numpy.array(extrem)
    exArrMax = exArr[1::2]
    maxVals = numpy.interp(exArrMax, T, fld)
    for i in range(exArrMax.shape[0]):
        pylab.semilogy(T[exArrMax[i]], fld[exArrMax[i]], 'ro')

    Tmax = T[exArrMax[:]]
    fldMax = fld[exArrMax[:]]
    return exArr, Tmax, fldMax

def calcOmegaGamma(T, fld):
    extrem, Tmax, fldMax = makePlot(T, fld)
    idxMin = extrem[0::2]

    # compute time-period of oscillations
    Tp = []
    for i in range(idxMin.shape[0]-1):
        Tp.append(T[idxMin[i+1]] - T[idxMin[i]])

    # compute least-square fit
    A = numpy.zeros((Tmax.shape[0], 2), numpy.float)
    A[:,0] = Tmax
    A[:,1] = 1
    fit = numpy.linalg.lstsq(A, numpy.log(fldMax))

    F = fit[0][0]*T + fit[0][1]
    V = pylab.exp(F)
    pylab.semilogy(T, V, 'k-')
    pylab.xlabel('Time')
    pylab.ylabel(r'Field Energy $\int |\nabla\phi |^2 dx$')
    pylab.savefig('s347-field-energy.png')
    return numpy.array(Tp), fit[0]

dat = pylab.loadtxt("fieldEnergy")
T = dat[:,0]
E = dat[:,1]
Tp, fit = calcOmegaGamma(T, E)
print Tp.sum()/Tp.shape[0], -fit[0]

#pylab.figure(2)
#pylab.plot(Tp)

pylab.show()
