from pylab import *
import math

d = loadtxt("err.txt")
dx = d[:,0]
err = d[:,1]

for i in range(1,dx.shape[0]):
    dxOrder = math.log(err[i]/err[i-1])/math.log(dx[i]/dx[i-1])
    print(dx[i], dxOrder)
