from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')

E0 = 0.5
w = 0.5

def wR(t):
    wx = E0/2*(t*cos(t)+sin(t))
    wy = -E0/2*t*sin(t)
    return wx, wy

Tex = linspace(0, 20, 500)
wx, wy = wR(Tex)

def plotFig():
    ux = zeros((21,), float)
    uy = zeros((21,), float)

    wx = zeros((21,), float)
    wy = zeros((21,), float)
    
    for i in range(21):
        print("Working on %d ..." % i)
        data = pg.GData("c3-oscc-E_ions_M1i_%d.bp" % i)
        dg = pg.data.GInterpModal(data, 2, "ms")
        XX, Ux = dg.interpolate(0)
        XX, Uy = dg.interpolate(1)

        ux[i] = Ux[0]
        uy[i] = Uy[0]        
        
        tm = data.time
        wx1, wy1 = wR(tm)
        wx[i] = wx1
        wy[i] = wy1

    return ux, uy

figure(1)
T = linspace(0, 20, 21)
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

savefig("c3-oscc-E-c-cmp.png", dpi=150)

show()    

