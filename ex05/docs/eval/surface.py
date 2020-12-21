import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3

plt.style.use("ggplot")


def oned(partitions, N):
    return (N - 1) * 2 * partitions


def twod(partitions, N, aspect):
    return ((np.sqrt(partitions) - 1) * 2 * (N * (1 + aspect)))


fig = plt.figure()

a = 1

pdata = np.arange(1, 16, 1)
ndata = np.arange(2, 32, 4)

ax = fig.add_subplot(111, projection="3d")
pdata, ndata = np.meshgrid(pdata, ndata)

odata = oned(pdata, ndata)
tdata = twod(pdata, ndata, a)

ax.plot_surface(pdata, ndata, odata)
ax.plot_surface(pdata, ndata, tdata)

ax.set_xlabel("Partitions")
ax.set_ylabel("Grid Dimension")
ax.set_zlabel("Number of Points to transfer")
plt.savefig(f"img/partition_{a}.pdf")
