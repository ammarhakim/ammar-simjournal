import os

for i in range(100):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s403-euler-implode-ds-2d_q_%d.h5 -c 0 --dont-show --save -o implode_%05d" % (i, i)
    print cmd
    os.system(cmd)
