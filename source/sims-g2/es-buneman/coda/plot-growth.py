from pylab import *
style.use('postgkyl.mplstyle')

M = linspace(20, 1836.2*45, 400)
g = sqrt(3.0)/2.0*(0.5/M)**(1/3)*(1-0.5*(0.5/M)**(1/3))

dat = loadtxt("growth-rates.txt")

figure(1)

plot(M, g, '-k', label=r'Max $\gamma$')
#semilogx(M, g, '-k', label=r'Max $\gamma$')
plot(dat[:,0], dat[:,1], 'ro', label=r'G2  $\gamma$')
#semilogx(dat[:,0], dat[:,1], 'ro', label=r'G2  $\gamma$')
grid()
xlabel(r'Mass Ratio $m_i/m_e$')
ylabel(r'Growth Rate $\gamma/\omega_{pe}$')
#title('Buneman Instability')
#legend(loc='best')

savefig('buneman-growth-cmp.png', dpi=150)
show()
