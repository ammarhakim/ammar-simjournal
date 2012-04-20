from pylab import *

d = loadtxt("../s122/totalEnstrophy")
plot(d[:,0], d[:,1], label='DG2-UF')
print "DG2-UF", (d[-1,1]-d[0,1])/d[0,1]*100

d = loadtxt("../s123/totalEnstrophy")
plot(d[:,0], d[:,1], label='DG2-CF')
print "DG2-CF", (d[-1,1]-d[0,1])/d[0,1]*100

d = loadtxt("../s124/totalEnstrophy")
plot(d[:,0], d[:,1], label='DG3-UF')
print "DG3-UF", (d[-1,1]-d[0,1])/d[0,1]*100

d = loadtxt("totalEnstrophy")
plot(d[:,0], d[:,1], label='DG2-CF')
print "DG3-CF", (d[-1,1]-d[0,1])/d[0,1]*100

legend()

title('Enstrophy history')
xlabel('Time [s]')
ylabel('Total Enstrophy')
savefig('enstrophy-hist-cmp.png')


