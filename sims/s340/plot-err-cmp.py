from pylab import *

solL2 = 3.1399614 # this is L2 norm of initial condition
baseList = ["../s337/s337-modal-dg-diffuse", "../s338/s338-modal-dg-diffuse", 
            "../s339/s339-modal-dg-diffuse", "../s340/s340-modal-dg-diffuse"]
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
savefig('s337-s338-s339-s340-rel-rms-err.png')

show()


