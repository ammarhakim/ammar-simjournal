import numpy as np
from pylab import *

style.use('postgkyl.mplstyle')

def buneman_k2w_cold(k, m):
    """
    Dispersion relation of two-fluid Bunenman instability in
    a cold plasma.
    
    Compute complex frequency w as roots of DR=0. DR is
    the dispersion relation equation rearranged as a
    polynommial equation where w is the unknown.
    The roots are ordered by imagninary parts (growth rate)
    first and then by real parts, both from large to small.
    
    Set se and si to zero for a cold plasma.
    
    Args:
        k: k*v0/wpe
        m: mi/me

    Returns:
        ws: Complex frequency for the given k.
    """
    ws = np.roots((m, -2 * k * m, k**2 * m - (m + 1), 2 * k,  -k**2))
    return ws[np.argsort(ws.imag + 1j * ws.real)[::-1]]

nk = 500
ks = np.linspace(0.01, 1.75, nk)

ws25 = np.zeros((nk, 4), dtype=np.complex128)
for i, k in enumerate(ks):
    ws25[i,:] = buneman_k2w_cold(k, 25.0)

ws200 = np.zeros((nk, 4), dtype=np.complex128)
for i, k in enumerate(ks):
    ws200[i,:] = buneman_k2w_cold(k, 200.0)

ws1836 = np.zeros((nk, 4), dtype=np.complex128)
for i, k in enumerate(ks):
    ws1836[i,:] = buneman_k2w_cold(k, 1836.2)

# plot imag part of dispersion relation
figure(1)
plot(ks, imag(ws25[:,0]), 'r.', markersize=2)
plot(ks, imag(ws25[:,1]), 'r.', markersize=2)
plot(ks, imag(ws25[:,2]), 'r.', markersize=2)
plot(ks, imag(ws25[:,3]), 'r.', markersize=2)

plot(ks, imag(ws200[:,0]), 'k.', markersize=2)
plot(ks, imag(ws200[:,1]), 'k.', markersize=2)
plot(ks, imag(ws200[:,2]), 'k.', markersize=2)
plot(ks, imag(ws200[:,3]), 'k.', markersize=2)

plot(ks, imag(ws1836[:,0]), 'm.', markersize=2)
plot(ks, imag(ws1836[:,1]), 'm.', markersize=2)
plot(ks, imag(ws1836[:,2]), 'm.', markersize=2)
plot(ks, imag(ws1836[:,3]), 'm.', markersize=2)

xlabel(r'$k V_0/\omega_{pe}$')
ylabel(r'$\gamma/\omega_{pe}$')
grid()
savefig('gamma-vs-k-cold.png', dpi=150)

figure(2)
plot(real(ws25[:,0]), imag(ws25[:,0]), 'r.', markersize=2)
plot(real(ws25[:,1]), imag(ws25[:,1]), 'r.', markersize=2)
plot(real(ws25[:,2]), imag(ws25[:,2]), 'r.', markersize=2)
plot(real(ws25[:,3]), imag(ws25[:,3]), 'r.', markersize=2)

plot(real(ws200[:,0]), imag(ws200[:,0]), 'k.', markersize=2)
plot(real(ws200[:,1]), imag(ws200[:,1]), 'k.', markersize=2)
plot(real(ws200[:,2]), imag(ws200[:,2]), 'k.', markersize=2)
plot(real(ws200[:,3]), imag(ws200[:,3]), 'k.', markersize=2)

plot(real(ws1836[:,0]), imag(ws1836[:,0]), 'm.', markersize=2)
plot(real(ws1836[:,1]), imag(ws1836[:,1]), 'm.', markersize=2)
plot(real(ws1836[:,2]), imag(ws1836[:,2]), 'm.', markersize=2)
plot(real(ws1836[:,3]), imag(ws1836[:,3]), 'm.', markersize=2)

axis('image')
gca().set_xlim([0.0, 0.5])
grid()
xlabel(r'$\omega/\omega_{pe}$')
ylabel(r'$\gamma/\omega_{pe}$')

savefig('omega-vs-gamma-cold.png', dpi=150)

show()

