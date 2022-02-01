import os

for i in range(100):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s301-5m-double-periodic_q_%d.h5 -c 3 --dont-show --save -o s301-5m-double-periodic_jz_%05d" % (i, i)
    print cmd
    os.system(cmd)
