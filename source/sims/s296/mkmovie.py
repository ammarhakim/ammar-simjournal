import os

for i in range(100):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s296-harris-tenmom_q_%d.h5 -c 3 --dont-show --save -o s296-harris-tenmom_jz_%05d" % (i, i)
    print cmd
    os.system(cmd)
