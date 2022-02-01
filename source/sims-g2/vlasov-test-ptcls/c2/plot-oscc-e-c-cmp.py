from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

E0 = 1.0
w = 0.5

def wNR(t):
    wx = E0/(1-w**2)*(sin(t)-w*sin(w*t))
    wy = E0/(1-w**2)*(cos(t)-cos(w*t))
    return wx, wy

Tex = linspace(0, 100, 5000)
wx, wy = wNR(Tex)

def plotFig():
    ux = zeros((101,), float)
    uy = zeros((101,), float)

    wx = zeros((101,), float)
    wy = zeros((101,), float)
    
    for i in range(101):
        print("Working on %d ..." % i)
        data = pg.GData("c2-oscc-E_ions_M1i_%d.bp" % i)
        dg = pg.data.GInterpModal(data, 2, "ms")
        XX, Ux = dg.interpolate(0)
        XX, Uy = dg.interpolate(1)

        ux[i] = Ux[0]
        uy[i] = Uy[0]        
        
        tm = data.time
        wx1, wy1 = wNR(tm)
        wx[i] = wx1
        wy[i] = wy1

    return ux, uy

figure(1)
T = linspace(0, 100, 101)
ux, uy = plotFig()
subplot(2,1,1)
plot(T, ux, 'ro', Tex, wx, 'k-')
ylabel('$u_x$')
grid()

subplot(2,1,2)
plot(T, uy, 'ro', Tex, wy, 'k-')
ylabel('$u_y$')
xlabel("Time [s]")
grid()

savefig("c2-oscc-E-c-cmp.png", dpi=150)

show()    

