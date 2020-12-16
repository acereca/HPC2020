import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as man

# figure and axes setup
fig, ax = plt.subplots()

harv = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                 [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                 [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                 [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                 [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                 [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
                 [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])

im = ax.imshow(harv, vmin=0, vmax=10)

# animation function, creates the plot frame by frame (i)
def animate(i):
    harv = np.random.random((7,7)) * 10
    im.set_data(harv)


anim = man.FuncAnimation(fig, animate, frames=120,
                               interval=1000/60, blit=False)
plt.tight_layout()
ax.set_xticks(np.arange(7))
ax.set_yticks(np.arange(7))
anim.save('img/heat.gif', fps=60, writer='imagemagick')
