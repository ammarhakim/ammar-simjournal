from pylab import *

solL2 = 3.1399614 # this is L2 norm of initial condition
baseList = ["../s341/s341-modal-dg-diffuse", "../s342/s342-modal-dg-diffuse", 
            "../s343/s343-modal-dg-diffuse", "../s344/s344-modal-dg-diffuse"]
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
savefig('s341-s342-s343-s344-rel-rms-err.png')

show()


