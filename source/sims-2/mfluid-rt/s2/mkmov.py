import os

for i in range(101):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s2-5m-2d-rt_q_%d.h5 -c 0 -t 'N_e t=%g' --dont-show --save -o s2-5m-2d-rt_numElc_%05d" % (i, 5*i, i)
    print cmd
    os.system(cmd)
