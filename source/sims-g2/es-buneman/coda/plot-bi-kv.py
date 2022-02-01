from pylab import *
style.use('postgkyl.mplstyle')

dat = loadtxt("bi-kv-1.txt")
plot(dat[:,0], dat[:,1], label=r'$kV_0/\omega_{pe}=1$')

dat = loadtxt("bi-kv-34.txt")
plot(dat[:,0], dat[:,1], label=r'$kV_0/\omega_{pe}=3/4$')

dat = loadtxt("bi-kv-12.txt")
plot(dat[:,0], dat[:,1], label=r'$kV_0/\omega_{pe}=1/2$')

dat = loadtxt("bi-kv-43.txt")
plot(dat[:,0], dat[:,1], label=r'$kV_0/\omega_{pe}=4/3$')

legend(loc='best')

xlabel(r'Mass Ratio $m_i/m_e$')
ylabel(r'Growth Rate $\gamma/\omega_{pe}$')

grid()

savefig('buneman-kv-cmp.png', dpi=150)

show()
