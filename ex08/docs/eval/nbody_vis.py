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


print(df[df['step'] == 0].iloc[0], df[df['step'] == df['step'].max()].iloc[0])

# fig, axs = plt.subplots(1, 3, figsize=(15,5))

# i = 0
# ims = []

# ax = axs[0]
# ims.append(
    # ax.errorbar(
        # df[df['step'] == i]['pos.x'],
        # df[df['step'] == i]['pos.y'],
        # fmt=".",
        # # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    # )[0]
# )
# ax.set_xlim([df['pos.x'].min(),df['pos.x'].max()])
# ax.set_ylim([df['pos.y'].min(),df['pos.y'].max()])

# ax = axs[1]
# ims.append(
    # ax.errorbar(
        # df[df['step'] == i]['pos.x'],
        # df[df['step'] == i]['pos.z'],
        # fmt=".",
        # # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    # )[0]
# )
# ax.set_xlim([df['pos.x'].min(),df['pos.x'].max()])
# ax.set_ylim([df['pos.z'].min(),df['pos.z'].max()])

# ax = axs[2]
# ims.append(
    # ax.errorbar(
        # df[df['step'] == i]['pos.y'],
        # df[df['step'] == i]['pos.z'],
        # fmt=".",
        # # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    # )[0]
# )
# ax.set_xlim([df['pos.y'].min(),df['pos.y'].max()])
# ax.set_ylim([df['pos.z'].min(),df['pos.z'].max()])

# # animation function, creates the plot frame by frame (i)
# def animate(i):
    # ims[0].set_data(df[df['step'] == i]['pos.x'], df[df['step'] == i]['pos.y'])
    # ims[1].set_data(df[df['step'] == i]['pos.x'], df[df['step'] == i]['pos.z'])
    # ims[2].set_data(df[df['step'] == i]['pos.y'], df[df['step'] == i]['pos.z'])
    # print(f"Progress: {round(100*i/df['step'].max())}% ({i}/{df['step'].max()})", end="\r")

# fps = 3
# anim = man.FuncAnimation(fig, animate, frames=df['step'].max(),
                               # interval=1000/fps, blit=False)
# plt.tight_layout()
# # ax.set_xticks([])
# # ax.set_yticks([])
# print("Stitching ...")
# anim.save('img/nbody_vis.mp4', fps=fps, dpi=300)
# print("Done!")

fig = plt.figure();
ax = fig.add_subplot(111, projection='3d')

i = 0
im, = ax.plot(
        df[df['step'] == i]['pos.x'],
        df[df['step'] == i]['pos.y'],
        df[df['step'] == i]['pos.z'],
        marker='o',
        linestyle='none'
        # ms=np.log(df[df['step'] == i]['mass [kg]']).tolist()
    )
# ax.set_xlim([df['pos.x'].min(),df['pos.x'].max()])
# ax.set_ylim([df['pos.y'].min(),df['pos.y'].max()])


# animation function, creates the plot frame by frame (i)
def animate(i):
    im.set_data(df[df['step'] == i]['pos.x'], df[df['step'] == i]['pos.y'])
    im.set_3d_properties(df[df['step'] == i]['pos.z'])
    print(f"Progress: {round(100*i/df['step'].max())}% ({i}/{df['step'].max()})", end="\r")

fps = df['step'].max()/20
anim = man.FuncAnimation(fig, animate, frames=df['step'].max(),
                               interval=1000/fps, blit=False)
plt.tight_layout()
# ax.set_xticks([])
# ax.set_yticks([])
print("Stitching ...")
anim.save('img/nbody_vis.mp4', fps=fps, dpi=300)
print("Done!")
