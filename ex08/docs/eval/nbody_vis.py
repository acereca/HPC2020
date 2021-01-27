import pandas as pd
import glob
import numpy as np
import argparse as ap
import matplotlib.pyplot as plt
from matplotlib import rc as rc
import matplotlib.animation as man

rc("text", usetex=True)
rc("font", **{"family": "serif", "serif": ["Computer Modern Roman"]})
plt.style.use("bmh")

df = pd.DataFrame()

for f in glob.glob("../code/out/bodies_r*"):
    df = df.append(
        pd.read_csv(
            f,
            sep=",",
            names=["step", "mass [kg]", 'pos.x', 'pos.y', 'pos.z'],
            comment="#"
        )
    )


print(df[df['step'] == 0].iloc[0], df[df['step'] == 199].iloc[0])
print(df[df['step'] == 0].iloc[1], df[df['step'] == 199].iloc[1])

fig, axs = plt.subplots(3, 1)

i = 0
ims = []

ax = axs[0]
ims.append(
    ax.errorbar(
        df[df['step'] == i]['pos.x'],
        df[df['step'] == i]['pos.y'],
        fmt=".",
        # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    )[0]
)

ax = axs[1]
ims.append(
    ax.errorbar(
        df[df['step'] == i]['pos.x'],
        df[df['step'] == i]['pos.z'],
        fmt=".",
        # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    )[0]
)

ax = axs[2]
ims.append(
    ax.errorbar(
        df[df['step'] == i]['pos.y'],
        df[df['step'] == i]['pos.z'],
        fmt=".",
        # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    )[0]
)

for ax in axs:
    ax.set_xlim([-1e2,1e2])
    ax.set_ylim([-1e2,1e2])
# animation function, creates the plot frame by frame (i)
def animate(i):
    ims[0].set_data(df[df['step'] == i]['pos.x'], df[df['step'] == i]['pos.y'])
    ims[1].set_data(df[df['step'] == i]['pos.x'], df[df['step'] == i]['pos.z'])
    ims[2].set_data(df[df['step'] == i]['pos.y'], df[df['step'] == i]['pos.z'])
    print(f"Progress: {round(100*i/df['step'].max())}% ({i}/{df['step'].max()})", end="\r")

fps = 3
anim = man.FuncAnimation(fig, animate, frames=df['step'].max(),
                               interval=1000/fps, blit=False)
plt.tight_layout()
# ax.set_xticks([])
# ax.set_yticks([])
print("Stitching ...")
anim.save('img/nbody_vis.mp4', fps=fps, dpi=300)
print("Done!")
