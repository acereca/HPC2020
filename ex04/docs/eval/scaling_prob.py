import pandas as pd
import glob
import argparse as ap
import numpy as np
import scipy.optimize as opt

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
            group['dim'],
            # group['GFLOPS64/s'],
            group['GFLOPS64/s'].mean(),
            # yerr=group['GFLOPS64/s'].std(),
            # fmt=".",
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
