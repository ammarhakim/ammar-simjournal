from pylab import *
import tables

# read data from various files
s122 = tables.openFile("../s122/s122-pb-advection-2d_totalEnstrophy_1.h5")
s123 = tables.openFile("../s123/s123-pb-advection-2d_totalEnstrophy_1.h5")
s124 = tables.openFile("s124-pb-advection-2d_totalEnstrophy_1.h5")

# make plots
plot(s122.root.DataStruct.timeMesh, s122.root.DataStruct.data, label='CFL 0.2')
plot(s123.root.DataStruct.timeMesh, s123.root.DataStruct.data, label='CFL 0.1')
plot(s124.root.DataStruct.timeMesh, s124.root.DataStruct.data, label='CFL 0.05')
xlabel('Time [s]')
ylabel('Enstrophy')
title('Total Enstrophy History, DG Order 2')
legend(loc = 3)
savefig('s124-rk2-enstrophy-history.png')

def tellError(dat):
    print (dat[-1]-dat[0])/dat[0]

tellError(s122.root.DataStruct.data)
tellError(s123.root.DataStruct.data)
tellError(s124.root.DataStruct.data)

show()
