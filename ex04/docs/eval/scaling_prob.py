import pandas as pd
import glob
import argparse as ap
import numpy as np
import scipy.optimize as opt
import VisTools.tex as vt

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpp
from matplotlib import rc as rc

rc("text", usetex=True)
rc("font", **{"family": "serif", "serif": ["Computer Modern Roman"]})
plt.style.use("ggplot")
# plt.figure(figsize=(8,6))

p = ap.ArgumentParser()
p.add_argument("slurm_job_id", metavar="J", type=int, nargs="?", default=83897)
args = p.parse_args()

df = pd.DataFrame()
df = df.append(
    pd.read_csv(
        "../code/out/"+str(args.slurm_job_id) + "",
        sep=",",
        names=["t/ns", "dim", "nodes", "GFLOPS64/s"],
        comment="#",
    )
)

# df = df.drop(df[df["s_msg/kB"].apply(type) != type(1)].index)

print(df)

with mpp.PdfPages("img/scaling_prob.pdf") as pdf:

    for name, group in df.groupby('nodes'):
        plt.errorbar(
            group.groupby('dim').groups.keys(),
            # group['GFLOPS64/s'],
            group.groupby('dim')['GFLOPS64/s'].mean(),
            yerr=group.groupby('dim')['GFLOPS64/s'].std(),
            fmt="--",
            marker=".",
            linewidth=1,
            elinewidth=.5,
            capsize=2,
            alpha=.8,
            label=name
        )

    plt.xlim([100,10000])
    plt.legend()
    plt.xscale("log")
    # plt.yscale("log")

    plt.title("Problem size dependent scaling of MMUL w/ MPI")
    plt.ylabel("GFLOPS64 / s")
    plt.xlabel("Matrix Dimension")

    plt.tight_layout()
    pdf.savefig()
    plt.close()

data = df.groupby(['dim', 'nodes']).agg({"GFLOPS64/s": ['mean', 'std']}).round(3)
data.to_latex(
    "data/scaling_prob.tex",
    escape=False,
    # formatters=[formatterfunc]*len(table.columns),
    # index=False,
    encoding='utf-8')

