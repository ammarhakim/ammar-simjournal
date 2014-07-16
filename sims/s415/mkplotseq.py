import os

for i in range(11):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s415-euler-nozzle-2d_q_%d.h5 -c 0 -m s415-euler-nozzle-2d_inOut.h5 -t 'Density t=%g' --dont-show --save -o s415-euler-nozzle-2d_%05d" % (i, i, i)
    print cmd
    os.system(cmd)
