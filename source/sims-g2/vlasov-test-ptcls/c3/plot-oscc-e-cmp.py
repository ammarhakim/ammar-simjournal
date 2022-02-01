from pylab import *
import postgkyl as pg

style.use('../code/postgkyl.mplstyle')
subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9,
                wspace=0.2, hspace=0.2)

E0 = 0.5

def wR(t):
    wx = E0/2*(t*cos(t)+sin(t))
    wy = -E0/2*t*sin(t)
    return wx, wy

T = linspace(0, 20, 500)
WX, WY = wR(T)

vx = linspace(-8, 8, 60)
vy = linspace(-8, 8, 60)
VX, VY = meshgrid(vx, vy)

def fexact(wx, wy):
    v2 = (VX-wx)**2 + (VY-wy)**2
    return 1.0/(2*pi)*exp(-v2/2)

def plotFig(nr,nc,i,fr, showX=False):
    print("Working on %d ..." % i)
    data = pg.GData("c3-oscc-E_ions_%d.bp" % fr)
    dg = pg.data.GInterpModal(data, 2, "ms")
    XX, q = dg.interpolate()

    tm = data.time
    wx, wy = wR(tm)

    qEx = fexact(wx, wy)

    f = subplot(nr,nc,i)    
    pcolormesh(XX[1], XX[2], transpose(q[3,:,:,0]))
    plot(WX, WY, 'w', linewidth=0.5)
    if showX:
        xlabel("$v_x$")
    
    ylabel('$v_y$')
    grid()        
    axis('image')

    f = subplot(nr,nc,i+1)
    pcolormesh(VX, VY, qEx)
    plot(WX, WY, 'w', linewidth=0.5)
    if showX:
        xlabel("$v_x$")
    grid()
    axis('image')


figure(1)    
plotFig(2,2,1, 10)
plotFig(2,2,3, 20, True)
savefig("c3-oscc-E-cmp.png", dpi=150)

show()    

