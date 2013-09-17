import os

for i in range(61):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s285-pulsebox-wave_q_%d.h5 -c 7 --dont-show --save -o s285-pulsebox-wave_ez_%05d" % (i, i)
    print cmd
    os.system(cmd)
