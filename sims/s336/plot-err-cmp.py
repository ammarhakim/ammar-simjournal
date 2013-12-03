from pylab import *

solL2 = 3.1399614 # this is L2 norm of initial condition
baseList = ["../s333/s333-modal-dg-diffuse", "../s334/s334-modal-dg-diffuse", 
            "../s335/s335-modal-dg-diffuse", "../s336/s336-modal-dg-diffuse"]
titleStr = ["LDG-L", "LDG-R", "LDG-S", "RDG"]

count = 0
figure(1)
for baseName in baseList:
    dat = loadtxt("%s_l2Norm.txt" % baseName)
    T = dat[:,0]
    val = dat[:,1]
    semilogy(T, sqrt(val/(exp(-T)*solL2)), label=titleStr[count])
    count = count+1
legend(loc='best')
title('Relative RMS Error')
savefig('s333-s334-s335-s336-rel-rms-err.png')

show()


