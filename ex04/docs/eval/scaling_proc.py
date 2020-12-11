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
p.add_argument("slurm_job_id", metavar="J", type=int, nargs="?", default=83896)
args = p.parse_args()

df = pd.DataFrame()
df = df.append(
    pd.read_csv(
        "../code/out/"+str(args.slurm_job_id),
        sep=",",
        names=["t/ns", "dim", "nodes", "GFLOPS64/s"],
        comment="#",
    )
)

# df = df.drop(df[df["s_msg/kB"].apply(type) != type(1)].index)

print(df)

with mpp.PdfPages("img/scaling_proc.pdf") as pdf:

    plt.errorbar(
        df.groupby('nodes').groups.keys(),
        df.groupby("nodes")['GFLOPS64/s'].mean(),
        # df['GFLOPS64/s'].mean(),
        yerr=df.groupby("nodes")['GFLOPS64/s'].std(),
        fmt="--",
        linewidth=1,
        marker="."
    )


    plt.legend()
    # plt.xscale("log")
    # plt.yscale("log")

    plt.title("Process dependent scaling of 2k MMUL w/ MPI")
    plt.ylabel("GFLOPS64 / s")
    plt.xlabel("Processes")

    plt.tight_layout()
    pdf.savefig()
    plt.close()

dfnodes = df.groupby('nodes')

data = dfnodes.agg({"GFLOPS64/s": ['mean', 'std']}).round(3)
data = pd.concat([data, (dfnodes.agg({'t/ns': ['mean', 'std']})/1e9).round(3)], axis=1)
data[('Speedup', "")] = (dfnodes['t/ns'].mean().iloc[0]/dfnodes['t/ns'].mean()).round(3)
data.to_latex(
    "data/scaling_proc.tex",
    escape=False,
    # columns=[("Time/s", ""), ("Speedup", ""), ("GFLOPS64/s", "mean"), ("GFLOPS64/s", "std")],
    # formatters=[(lambda t: t/1e9), (lambda f: round(f, 3)), None],
    multicolumn=True,
    encoding='utf-8')

