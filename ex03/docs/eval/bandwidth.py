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
p.add_argument("slurm_job_id1", metavar="J", type=int, nargs="?", default=82834)
p.add_argument("slurm_job_id2", metavar="J", type=int, nargs="?", default=82835)
args = p.parse_args()

df = pd.DataFrame()
for f in glob.glob("../code/out/" + str(args.slurm_job_id1) + "_*"):
    df = df.append(
        pd.read_csv(
            f,
            sep=",",
            names=["t_msg/ns", "s_msg/kB", "nodes", "blocking"],
            # dtype={
            # "t_msg/ns": np.int64,
            # "s_msg/kB": np.int32,
            # "nodes": np.int32,
            # "blocking": np.float64,
            # },
            comment="#",
        )
    )

for f in glob.glob("../code/out/" + str(args.slurm_job_id2) + "_*"):
    df = df.append(
        pd.read_csv(
            f,
            sep=",",
            names=["t_msg/ns", "s_msg/kB", "nodes", "blocking"],
            dtype={
                # "t_msg/ns": np.int64,
                # "s_msg/kB": np.int32,
                # "nodes": np.int32,
                # "blocking": np.float64,
            },
            comment="#",
        )
    )

df = df.drop(df[df["s_msg/kB"].apply(type) != type(1)].index)

# df = df.drop(df[df['t_msg/ns'] > 1e7].index)
print(df)

with mpp.PdfPages("img/bandwidth.pdf") as pdf:

    for (title, block), group in df.groupby(["nodes", "blocking"]):
        print((title, block))
        sizes = group["s_msg/kB"].unique()
        t_means = np.array(
            [group[group["s_msg/kB"] == s]["t_msg/ns"].mean() for s in sizes]
        )
        means = (sizes * 1024) / (t_means / 1e9)
        t_stdds = np.array(
            [group[group["s_msg/kB"] == s]["t_msg/ns"].std() for s in sizes]
        )
        stdds = (sizes * 1024) / (t_stdds / 1e9)
        plt.errorbar(
            sizes,
            means / 1024 / 1024,
            # yerr=stdds,
            fmt=".",
            alpha=0.8,
            ms=6,
            elinewidth=1,
            label=f"{title} Node{'s' if title == 2 else ''}{', blocking send' if block else ''}",
        )

        m = np.max(means)/1024/1024

        plt.axhline(
            m,
            linewidth=1,
            label=r"${:.1f}$ MB/s".format(m),
            c=f"C{(title-1)*2+block}",
            ls="dashed",
        )
    plt.legend()
    # plt.xscale("log")
    plt.yscale("log")

    plt.title("Bandwidth for different Message Sizes")
    plt.ylabel("Bandwidth / MB/s")
    plt.xlabel("Message Size / kB")

    plt.tight_layout()
    pdf.savefig()
    plt.close()
