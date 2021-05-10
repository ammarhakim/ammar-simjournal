import postgkyl as pg
from pylab import *

count = 1

data = pg.GData("s2-kh-euler-euler_20.gkyl")
s2, iso_s2  = pg.diagnostics.fft(data, psd=True, iso=True)
rho_spec = iso_s2[:,0]

for i in range(21,41):
    data = pg.GData("s2-kh-euler-euler_%d.gkyl" % i)
    s2, iso_s2 = pg.diagnostics.fft(data, psd=True, iso=True)
    rho_spec = rho_spec + iso_s2[:,0]
    count = count + 1

loglog(s2[0], rho_spec/count)
show()
