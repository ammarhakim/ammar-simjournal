import os

for i in range(70):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s1-5m-2d-rt_q_%d.h5 -c 0 -t 'Electron density' --dont-show --save -o s1-5m-2d-rt_numElc_%05d" % (i, i)
    print cmd
    os.system(cmd)
