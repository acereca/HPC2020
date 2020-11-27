import pandas as pd
import glob
import argparse as ap
import numpy as np
import scipy.optimize as opt

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpp
from matplotlib import rc as rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
plt.style.use('bmh')
# plt.figure(figsize=())

p = ap.ArgumentParser()
p.add_argument("slurm_job_id1", metavar="J", type=int, nargs="?")
p.add_argument("slurm_job_id2", metavar="J", type=int, nargs="?")
args = p.parse_args()

df = pd.DataFrame()
for f in glob.glob("../code/out/" + str(args.slurm_job_id1) + "_*"):
    df = df.append(pd.read_csv(f, sep=",", names=["t_msg/ns", "s_msg/kB", "nodes"], comment="#"))

for f in glob.glob("../code/out/" + str(args.slurm_job_id2) + "_*"):
    df = df.append(pd.read_csv(f, sep=",", names=["t_msg/ns", "s_msg/kB", "nodes"], comment="#"))

print(df)
# df = df.drop(df[df['t_msg/ns'] > 1e7].index)
print(df)

with mpp.PdfPages('img/latency.pdf') as pdf:

    for title, group in df.groupby('nodes'):
        sizes = group['s_msg/kB'].unique()
        means = [group[group['s_msg/kB'] == s]['t_msg/ns'].mean()/1e3 for s in sizes]
        stdds = [group[group['s_msg/kB'] == s]['t_msg/ns'].std()/1e3 for s in sizes]
        plt.errorbar(
            sizes,
            means,
            yerr=stdds,
            fmt=".",
            alpha=.8,
            ms=5,
            elinewidth=1,
            label=f"{title} Node{'s' if title == 2 else ''}"
        )

        popt, pcov = opt.curve_fit(
            lambda x, m, c: x*m+c,
            sizes,
            means,
            p0=[10, 0],
            sigma=stdds,
            absolute_sigma=True
        )

        plt.plot(
            sorted(sizes),
            popt[0]*np.array(sorted(sizes))+popt[1],
            linewidth=1,
            label=r"${:.1f} s_{{msg}} + {:.1f}$".format(*popt),
            c=f"C{title-1}",
            ls="dashed"
        )
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')

    plt.title("Roundtrip Times for different Message Sizes")
    plt.ylabel("Roundtrip Time / ms")
    plt.xlabel("Message Size / kB")

    pdf.savefig()
    plt.close()
