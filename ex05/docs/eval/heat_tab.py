import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("ggplot")

tpit = "Time per iteration/s"
fpit = "Flops per iteration"

data = pd.read_csv("../code/out/85027", names=["dim", "iterations", tpit])

data[tpit] /= 1e9
data[fpit] = 7 * (data['dim'] - 1) * (data['dim'] - 1)
data['Flops total'] = data['iterations'] * data[fpit]
data['GFLOP64/s'] = data[fpit] / data[tpit] / 1e9
data['Grid Size'] = data['dim'].astype(str) + "x" + data['dim'].astype(str)

data_out = data[['Grid Size', tpit, 'Flops total', 'GFLOP64/s']]

data_out.to_latex(
    "data/heat.tex",
    float_format="%.4f",
    index=False,
    escape=False,
    encoding="utf-8"
)

plt.errorbar(
    data['dim'],
    data['GFLOP64/s'],
    fmt="--",
    marker=".",
    # ms=2
    label="GFLOPS64"
)

plt.axvline(np.sqrt(10e6/8)/2, ls="dotted", label="Max to fit in Cache")

plt.xscale("log")
plt.xlabel("Grid dimension")
plt.ylabel("GFLOPS64/s")
plt.tight_layout()
plt.legend()
plt.savefig("img/heat.pdf")
