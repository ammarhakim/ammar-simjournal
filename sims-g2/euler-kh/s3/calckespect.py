from pylab import *
import postgkyl

style.use('postgkyl.mplstyle')

d = postgkyl.GData("s3-euler-kh_ke-spect_00190.bp")
s = d.getValues()
g = d.getGrid()
dx = g[0][1]-g[0][0]
gc = linspace(g[0][0]+0.5*dx, g[0][-1]-0.5*dx, g[0].shape[0]-1)

for i in range(191, 201):
    d = postgkyl.GData("s3-euler-kh_ke-spect_00%d.bp" % i)
    s = s + d.getValues()

s = s/11.0
loglog(gc, s)
grid()

x1 = 10.632
y1 = 32.86
x2 = 105.244
y2 = 0.011815
x3 = 990.384
y3 = 6.11e-9

p1 = log(y2/y1)/log(x2/x1)
p2 = log(y3/y2)/log(x3/x2)
print(p1, p2)

X1 = linspace(x1, x2, 20)
Y1 = 100*(X1/X1[0])**p1
loglog(X1, Y1, 'r-')
text(50,10,r'$k^{-3.5}$')

X2 = linspace(x2, x3, 20)
Y2 = 0.75e-1*(X2/X2[0])**p2
loglog(X2, Y2, 'r-')
text(500,1e-4,r'$k^{-6.5}$')

xlabel(r'$k/2\pi$')
ylabel(r'$KE$')
savefig('ke-spect.png')

show()
    
