import os

for i in range(21):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p gemguide-5m_q_%d.h5 -c 3 --dont-show --save -o gemguide-5m_%05d" % (i, i)
    print cmd
    os.system(cmd)
