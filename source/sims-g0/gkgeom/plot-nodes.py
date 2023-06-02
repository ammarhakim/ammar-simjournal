import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import postgkyl as pg

data = pg.GData("gkgeom_node_coords.gkyl")
vals = data.getValues()
R = vals[:,:,0]
Z = vals[:,:,1]


#plt.plot(R,Z,marker=".", color="k", linestyle="none")
plt.scatter(R,Z, marker=".")
segs1 = np.stack((R,Z), axis=2)
segs2 = segs1.transpose(1,0,2)
plt.gca().add_collection(LineCollection(segs1))
plt.gca().add_collection(LineCollection(segs2))
plt.axis("equal")
plt.show()
