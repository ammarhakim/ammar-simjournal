import exceptions
import numpy

class ExtractFluidVars1D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'xx' : 4+f, 'xy' : 5+f, 'xz' : 6+f,
                  'yy' : 7+f, 'yz' : 8+f,
                  'zz' : 9+f}

    def getRho(self, q):
        return q[:,self.v['n']]

    def getU(self, q):
        return q[:,self.v['u']]/self.getRho(q)

    def getV(self, q):
        return q[:,self.v['v']]/self.getRho(q)

    def getW(self, q):
        return q[:,self.v['w']]/self.getRho(q)

    def getPxx(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        return q[:,self.v['xx']] - r*u*u

    def getPxy(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        v = self.getV(q)
        return q[:,self.v['xy']] - r*u*v

    def getPxz(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        w = self.getW(q)
        return q[:,self.v['xz']] - r*u*w

    def getPyy(self, q):
        r = self.getRho(q)
        v = self.getV(q)
        return q[:,self.v['yy']] - r*v*v

    def getPyz(self, q):
        r = self.getRho(q)
        v = self.getV(q)
        w = self.getW(q)
        return q[:,self.v['yz']] - r*v*w

    def getPzz(self, q):
        r = self.getRho(q)
        w = self.getW(q)
        return q[:,self.v['zz']] - r*w*w

    def getTemp(self, q):
        return 1/3.0*(self.getPxx(q)+self.getPyy(q)+self.getPzz(q))/(self.getRho(q)/mass)

    def getThermVel(self, q):
        return numpy.sqrt(self.getTemp(q)/mass)

class ExtractEmVars1D(object):
    def __init__(self):
        pass

    def getEx(self, q):
        return q[:,20]
    def getEy(self, q):
        return q[:,21]
    def getEz(self, q):
        return q[:,22]    
    def getBx(self, q):
        return q[:,23]
    def getBy(self, q):
        return q[:,24]
    def getBz(self, q):
        return q[:,25]

class ExtractFluidVars2D(object):
    def __init__(self, offset):
        f = offset
        self.v = {'n' : 0+f, 'u' : 1+f, 'v' : 2+f, 'w' : 3+f,
                  'xx' : 4+f, 'xy' : 5+f, 'xz' : 6+f,
                  'yy' : 7+f, 'yz' : 8+f,
                  'zz' : 9+f}

    def getRho(self, q):
        return q[:,:,self.v['n']]

    def getU(self, q):
        return q[:,:,self.v['u']]/self.getRho(q)

    def getV(self, q):
        return q[:,:,self.v['v']]/self.getRho(q)

    def getW(self, q):
        return q[:,:,self.v['w']]/self.getRho(q)

    def getPxx(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        return q[:,:,self.v['xx']] - r*u*u

    def getPxy(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        v = self.getV(q)
        return q[:,:,self.v['xy']] - r*u*v

    def getPxz(self, q):
        r = self.getRho(q)
        u = self.getU(q)
        w = self.getW(q)
        return q[:,:,self.v['xz']] - r*u*w

    def getPyy(self, q):
        r = self.getRho(q)
        v = self.getV(q)
        return q[:,:,self.v['yy']] - r*v*v

    def getPyz(self, q):
        r = self.getRho(q)
        v = self.getV(q)
        w = self.getW(q)
        return q[:,:,self.v['yz']] - r*v*w

    def getPzz(self, q):
        r = self.getRho(q)
        w = self.getW(q)
        return q[:,:,self.v['zz']] - r*w*w

    def getTemp(self, q):
        return 1/3.0*(self.getPxx(q)+self.getPyy(q)+self.getPzz(q))/(self.getRho(q)/mass)

    def getThermVel(self, q):
        return numpy.sqrt(self.getTemp(q)/mass)

class ExtractEmVars2D(object):
    def getEx(self, q):
        return q[:,:,20]
    def getEy(self, q):
        return q[:,:,21]
    def getEz(self, q):
        return q[:,:,22]    
    def getBx(self, q):
        return q[:,:,23]
    def getBy(self, q):
        return q[:,:,24]
    def getBz(self, q):
        return q[:,:,25]

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
    def getPxx(self, q):
        return self.get(q, 'getPxx')
    def getPxy(self, q):
        return self.get(q, 'getPxy')
    def getPxz(self, q):
        return self.get(q, 'getPxz')
    def getPyy(self, q):
        return self.get(q, 'getPyy')
    def getPyz(self, q):
        return self.get(q, 'getPyz')
    def getPzz(self, q):
        return self.get(q, 'getPzz')
    def getTemp(self, q):
        return self.get(q, 'getTemp')
    def getThermVel(self, q):
        return self.get(q, 'getThermVel')

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
ionEx = ExtractFluidVars(10)
emEx = ExtractEmVars()

transformRegistry = {
    'rhoElc' : elcEx.getRho,
    'uElc' : elcEx.getU,
    'vElc' : elcEx.getV,
    'wElc' : elcEx.getW,
    'pxxElc' : elcEx.getPxx,
    'pxyElc' : elcEx.getPxy,
    'pxzElc' : elcEx.getPxz,
    'pyyElc' : elcEx.getPyy,
    'pyzElc' : elcEx.getPyz,
    'pzzElc' : elcEx.getPzz,
    'tempElc' : elcEx.getTemp,
    'thermalVelElc' : elcEx.getThermVel,
    'rhoIon' : ionEx.getRho,
    'uIon' : ionEx.getU,
    'vIon' : ionEx.getV,
    'wIon' : ionEx.getW,
    'pxxIon' : ionEx.getPxx,
    'pxyIon' : ionEx.getPxy,
    'pxzIon' : ionEx.getPxz,
    'pyyIon' : ionEx.getPyy,
    'pyzIon' : ionEx.getPyz,
    'pzzIon' : ionEx.getPzz,
    'tempIon' : ionEx.getTemp,
    'thermalVelIon' : ionEx.getThermVel,
    'Ex' : emEx.getEx,
    'Ey' : emEx.getEy,
    'Ez' : emEx.getEz,
    'Bx' : emEx.getBx,
    'By' : emEx.getBy,
    'Bz' : emEx.getBz,
}
