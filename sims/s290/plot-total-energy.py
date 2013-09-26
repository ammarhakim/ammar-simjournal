from pylab import *
import numpy
import tables

Lx = 50.0
Ly = 25.0
nFrame = 60
elcEnergy = numpy.zeros((nFrame+1,), numpy.float)
ionEnergy = numpy.zeros((nFrame+1,), numpy.float)
emEnergy = numpy.zeros((nFrame+1,), numpy.float)

Tm = linspace(0, nFrame, nFrame+1)

for i in range(0,nFrame+1):
    print "Working on %d .." % i
    fh = tables.openFile("s290-harris-tenmom_q_%d.h5" % i)
    q = fh.root.StructGridField

    dx = Lx/q.shape[0]
    dy = Ly/q.shape[1]

    elcEnergy[i] = 0.5*dx*dy*sum(q[:,:,4]+q[:,:,7]+q[:,:,9])
    ionEnergy[i] = 0.5*dx*dy*sum(q[:,:,14]+q[:,:,17]+q[:,:,19])
    emEnergy[i] = 0.5*dx*dy*sum(q[:,:,20]**2 + q[:,:,21]**2 + q[:,:,22]**2 + q[:,:,23]**2 + q[:,:,24]**2 + q[:,:,25]**2)

figure(1)
totalEnergy = elcEnergy+ionEnergy+emEnergy
plot(Tm, elcEnergy+ionEnergy+emEnergy, '-k', label='Total')
print "Energy lost is ", (totalEnergy[-1]-totalEnergy[0])/totalEnergy[0]*100, " percent"

figure(2)
subplot(2,2,1)
plot(Tm, elcEnergy, '-k', label='Elc')
legend(loc='best')
subplot(2,2,2)
plot(Tm, ionEnergy, '-r', label='Ion')
legend(loc='best')
subplot(2,2,3)
plot(Tm, emEnergy, '-b', label='EM')
legend(loc='best')

show()


