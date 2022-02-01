import exceptions

class ExtractFluidVars1D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'p' : 4+f}
        self.g = 5.0/3.0 # this should somehow be passed to this class

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

class ExtractEmVars1D(object):
    def __init__(self):
        pass

    def getEx(self, q):
        return q[:,10]
    def getEy(self, q):
        return q[:,11]
    def getEz(self, q):
        return q[:,12]    
    def getBx(self, q):
        return q[:,13]
    def getBy(self, q):
        return q[:,14]
    def getBz(self, q):
        return q[:,15]

class ExtractFluidVars2D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'p' : 4+f}
        self.g = 5.0/3.0 # this should somehow be passed to this class

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

class ExtractEmVars2D(object):
    def __init__(self):
        pass

    def getEx(self, q):
        return q[:,:,10]
    def getEy(self, q):
        return q[:,:,11]
    def getEz(self, q):
        return q[:,:,12]    
    def getBx(self, q):
        return q[:,:,13]
    def getBy(self, q):
        return q[:,:,14]
    def getBz(self, q):
        return q[:,:,15]    

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
    def getT(self, q):
        return self.getP(q)/(self.getRho(q)/mass)

class ExtractEmVars(object):
    def __init__(self):
        self.ex1D = ExtractEmVars1D()
        self.ex2D = ExtractEmVars2D()

    def getDim(self, q):
        return len(q.shape)-1

    def get(self, q, nm):
        if self.getDim(q) == 1:
            return ExtractEmVars1D.__dict__[nm](self.ex1D, q)
        elif self.getDim(q) == 2:
            return ExtractEmVars2D.__dict__[nm](self.ex2D, q)
        else:
            raise exceptions.RuntimeError("3D is not yet supported!")        

    def getEx(self, q):
        return self.get(q, 'getEx')
    def getEy(self, q):
        return self.get(q, 'getEy')
    def getEz(self, q):
        return self.get(q, 'getEz')
    def getBx(self, q):
        return self.get(q, 'getBx')
    def getBy(self, q):
        return self.get(q, 'getBy')
    def getBz(self, q):
        return self.get(q, 'getBz')

elcEx = ExtractFluidVars(0)
ionEx = ExtractFluidVars(5)
emEx = ExtractEmVars()

transformRegistry = {
    'rhoElc' : elcEx.getRho,
    'uElc' : elcEx.getU,
    'vElc' : elcEx.getV,
    'wElc' : elcEx.getW,
    'pElc' : elcEx.getP,
    'tempElc' : elcEx.getT,
    'rhoIon' : ionEx.getRho,
    'uIon' : ionEx.getU,
    'vIon' : ionEx.getV,
    'wIon' : ionEx.getW,
    'pIon' : ionEx.getP,
    'tempIon' : ionEx.getT,
    'Ex' : emEx.getEx,
    'Ey' : emEx.getEy,
    'Ez' : emEx.getEz,
    'Bx' : emEx.getBx,
    'By' : emEx.getBy,
    'Bz' : emEx.getBz,
}
