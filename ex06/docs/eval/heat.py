import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as man

# figure and axes setup
fig, ax = plt.subplots()

print("Loading ...")
data = pd.read_csv("../code/out/gall.csv", sep=";", comment="#", names=["iteration", "dimx", "dimy", "data"])
print("Done!")
data['data'] = data['data'].apply(lambda x: np.array([float(e) for e in x[2:-1].split(",")]))
# print(data['data'].iloc[0])

im = ax.imshow(np.resize(data.iloc[0]['data'], (data.iloc[0]['dimx'], data.iloc[0]['dimy'])), vmin=0, vmax=127)

# animation function, creates the plot frame by frame (i)
def animate(i):
    harv = np.resize(data.iloc[i]['data'], (data.iloc[i]['dimx'], data.iloc[i]['dimy']))
    im.set_data(harv)
    print(f"Progress: {round(100*i/data['iteration'].max())}% ({i}/{data['iteration'].max()})", end="\r")

fps = 30
anim = man.FuncAnimation(fig, animate, frames=data['iteration'].max(),
                               interval=1000/fps, blit=False)
plt.tight_layout()
ax.set_xticks([])
ax.set_yticks([])
print("Stitching ...")
anim.save('img/heat.mp4', fps=fps)
print("Done!")
