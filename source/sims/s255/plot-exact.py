import tables
from pylab import *
from scipy import special
erf = special.erf

X = linspace(0, 10, 64)
phi = tanh(10-X)
e = 1.0

def numDens(phi):
    return (exp(e*phi)*erf(sqrt(e*phi))-exp(e*phi))/2.0

fh = tables.openFile("s255-vlasov-fp-bounded_numDensity_10.h5")
q = fh.root.StructGridField

nx = numDens(phi)
figure(1)
plot(-phi, nx, '-r')
figure(2)
plot(X, nx, '-r', X, q[:,0], '-k')
axis('tight')
show()
