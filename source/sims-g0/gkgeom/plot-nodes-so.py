import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import postgkyl as pg

data = pg.GData("solovev_node_coords.gkyl")
vals = data.getValues()
R = vals[:,:,0]
Z = vals[:,:,1]

psid = pg.GData("solovev_psi.gkyl")
interp = pg.GInterpModal(psid,2,"ms")
grid, psi = interp.interpolate()
#pg.output.plot(psid, contour=True)

for d in range(len(grid)):
    grid[d] = 0.5*(grid[d][:-1] + grid[d][1:])

clevels = np.linspace(0.0496719, 0.1, 17)
plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="k")

clevels = np.linspace(0.0, 0.0497, 5)
plt.contour(grid[0], grid[1], psi[:,:,0].transpose(), levels=clevels, colors="r")

plt.plot(R,Z,marker=".", color="k", linestyle="none")
plt.scatter(R,Z, marker=".")
segs1 = np.stack((R,Z), axis=2)
segs2 = segs1.transpose(1,0,2)
plt.gca().add_collection(LineCollection(segs1))
plt.gca().add_collection(LineCollection(segs2))
plt.grid()
plt.axis("tight")
plt.axis("image")
plt.show()
