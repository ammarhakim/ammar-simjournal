from pylab import *

X = linspace(-1, 1, 50)
XX, YY = meshgrid(X, X)

def calcf(XX, YY, mu1):
    f1 = 0.5
    f2 = 1/(2*sqrt(3)*mu1)
    f3 = 1/(2*sqrt(3)*mu1)
    f4 = 1/(6*mu1**2)

    return f1*0.5 + f2*sqrt(3)*XX/2 + f3*sqrt(3)*YY/2 + 3*XX*YY/2

mu1 = 3.0/5.0
f1 = calcf(XX, YY, mu1)
pcolormesh(XX, YY, transpose(f1))
axis('image')
colorbar()
print("Max: %g. Min = %g" % (f1.max(), f1.min()))
show()
