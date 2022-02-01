import matplotlib
matplotlib.use('Agg')
from pylab import *
import tables
import vistools
import math
import tftenmoment

gasGamma = 5./3.
elcCharge = -1.0
ionCharge = 1.0
ionMass = 1.0
elcMass = ionMass/100.0
epsilon0 = 1.0
mu0 = 1.0
lightSpeed = 1/math.sqrt(epsilon0*mu0)

def getDx(grid):
    xl, yl = grid._v_attrs.vsLowerBounds
    xu, yu = grid._v_attrs.vsUpperBounds
    nx, ny = grid._v_attrs.vsNumCells
    dx = (xu-xl)/nx
    dy = (yu-yl)/ny
    return dx, dy

def getFluidEnergy(e, q):
    rho = e.getRho(q)
    u = e.getU(q)
    v = e.getV(q)
    w = e.getW(q)
    Pxx = e.getPxx(q)
    Pyy = e.getPyy(q)
    Pzz = e.getPzz(q)

    ke = 0.5*rho*(u*u+v*v+w*w)
    ie = 1.5*(Pxx+Pyy+Pzz)/3.0
    return ke, ie

def getEmEnergy(e, q):
    Ex = e.getEx(q)
    Ey = e.getEy(q)
    Ez = e.getEz(q)
    Bx = e.getBx(q)
    By = e.getBy(q)
    Bz = e.getBz(q)

    return 0.5*(Ex*Ex+Ey*Ey+Ez*Ez) + 0.5*(Bx*Bx+By*By+Bz*Bz)

def getBxByEnergy(e, q):
    Bx = e.getBx(q)
    By = e.getBy(q)

    return 0.5*(Bx*Bx+By*By)

def writeDat(fName, field):
    fp = open(fName, "w")
    savetxt(fp, field)
    fp.close()

totKe_e = []
totIe_e = []
totKe_i = []
totIe_i = []
totEm = []
totBxBy = []

for i in range(0, 101):
    print "Working on %d ..." % i
    fh = tables.openFile("s9-10m-karim-kh_q_%d.h5" % i)
    dx, dy = getDx(fh.root.StructGrid)
    q = fh.root.StructGridField

    # electrons
    ke, ie  = getFluidEnergy(tftenmoment.elcEx, q)
    totKe_e.append( dx*dy*sum(ke) )
    totIe_e.append( dx*dy*sum(ie) )

    # ions
    ke, ie  = getFluidEnergy(tftenmoment.ionEx, q)
    totKe_i.append( dx*dy*sum(ke) )
    totIe_i.append( dx*dy*sum(ie) )

    # EM
    em  = getEmEnergy(tftenmoment.emEx, q)
    totEm.append( dx*dy*sum(em) )

    # inplane B
    bxby  = getBxByEnergy(tftenmoment.emEx, q)
    totBxBy.append( dx*dy*sum(bxby) )

    fh.close()

# write data to file
writeDat("s9-10m-karim-kh_totKe_e.txt", totKe_e)
writeDat("s9-10m-karim-kh_totIe_e.txt", totIe_e)
writeDat("s9-10m-karim-kh_totKe_i.txt", totKe_i)
writeDat("s9-10m-karim-kh_totIe_i.txt", totIe_i)
writeDat("s9-10m-karim-kh_totEm.txt", totEm)
writeDat("s9-10m-karim-kh_totBxBy.txt", totBxBy)
