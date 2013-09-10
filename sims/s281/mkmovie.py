import os

for i in range(61):
    cmd = "/Users/ahakim/research/gkeyll-project/gkeyllall/gkeyll/scripts/gkeplot.py -p s281-gem-tenmom_q_%d.h5 -c 3 --dont-show --save -o s281-gem-tenmom_jz_%05d" % (i, i)
    print cmd
    os.system(cmd)
