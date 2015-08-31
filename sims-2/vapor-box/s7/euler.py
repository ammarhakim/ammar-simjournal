import exceptions
import numpy

gasGamma = 5.0/3.0

class ExtractFluidVars1D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'p' : 4+f}
        self.g = gasGamma # this should somehow be passed to this class

    def getRho(self, q):
        return q[:,self.v['n']]

    def getU(self, q):
        return q[:,self.v['u']]/self.getRho(q)

    def getV(self, q):
        return q[:,self.v['v']]/self.getRho(q)

    def getW(self, q):
        return q[:,self.v['w']]/self.getRho(q)

    def getP(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        v = self.getV(q)
        w = self.getW(q)
        return (q[:,self.v['p']] - 0.5*r*(u*u+v*v+w*w))*(self.g-1)

    def getMach(self, q):
        cs = numpy.sqrt(self.g*self.getP(q)/self.getRho(q))
        uu = numpy.sqrt(self.getU(q)**2 + self.getV(q)**2 + self.getW(q)**2)
        return uu/cs

    def getEnthalphy(self, q):
        pr = getP(q)
        return (q[:,self.v['p']])

    def getIe(self, q):
        r = self.getRho(q)
        return self.getP(q)/r/(self.g-1)

class ExtractFluidVars2D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'p' : 4+f}
        self.g = gasGamma # this should somehow be passed to this class

    def getRho(self, q):
        return q[:,:,self.v['n']]

    def getU(self, q):
        return q[:,:,self.v['u']]/self.getRho(q)

    def getV(self, q):
        return q[:,:,self.v['v']]/self.getRho(q)

    def getW(self, q):
        return q[:,:,self.v['w']]/self.getRho(q)

    def getP(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        v = self.getV(q)
        w = self.getW(q)
        return (q[:,:,self.v['p']] - 0.5*r*(u*u+v*v+w*w))*(self.g-1)

    def getMach(self, q):
        cs = numpy.sqrt(self.g*self.getP(q)/self.getRho(q))
        uu = numpy.sqrt(self.getU(q)**2 + self.getV(q)**2 + self.getW(q)**2)
        return uu/cs

    def getIe(self, q):
        r = self.getRho(q)
        return self.getP(q)/r/(self.g-1)

class ExtractFluidVars(object):
    def __init__(self, offset):
        self.ex1D = ExtractFluidVars1D(offset)
        self.ex2D = ExtractFluidVars2D(offset)

    def getDim(self, q):
        return len(q.shape)-1

    def get(self, q, nm):
        if self.getDim(q) == 1:
            return ExtractFluidVars1D.__dict__[nm](self.ex1D, q)
        elif self.getDim(q) == 2:
            return ExtractFluidVars2D.__dict__[nm](self.ex2D, q)
        else:
            raise exceptions.RuntimeError("3D is not yet supported!")

    def getRho(self, q):
        return self.get(q, 'getRho')
    def getU(self, q):
        return self.get(q, 'getU')
    def getV(self, q):
        return self.get(q, 'getV')
    def getW(self, q):
        return self.get(q, 'getW')
    def getP(self, q):
        return self.get(q, 'getP')
    def getMach(self, q):
        return self.get(q, 'getMach')
    def getIe(self, q):
        return self.get(q, 'getIe')

fluidEx = ExtractFluidVars(0)

transformRegistry = {
    'rho' : fluidEx.getRho,
    'u' : fluidEx.getU,
    'v' : fluidEx.getV,
    'w' : fluidEx.getW,
    'p' : fluidEx.getP,
    'ie' : fluidEx.getIe,
    'mach' : fluidEx.getMach
}
