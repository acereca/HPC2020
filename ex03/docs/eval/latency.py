import pandas as pd
import glob
import argparse as ap
import matplotlib.pyplot as plt
from matplotlib import rc as rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
plt.style.use('bmh')
# plt.figure(figsize=())

p = ap.ArgumentParser()
p.add_argument("slurm_job_id", metavar="J", type=int, nargs="?")
args = p.parse_args()

if args.slurm_job_id[0]

df = pd.DataFrame()
for f in glob.glob("out/" + str(args.slurm_job_id[0]) + "_*"):
    df = df.append(pd.read_csv(f, sep=",", names=["t_msg/ns", "s_msg/kB", "nodes"], comment="#"))

print(df)

plt.hist(df['t_msg/ns'], 100)


# for numm in df['nummsg'].unique():
    # f = (df['nummsg'] == numm)

    # plt.errorbar(
        # df['nproc'].unique(),
        # [df[f & (df['nproc'] == s)]['t_msg/ns'].mean()/1000 for s in df['nproc'].unique()],
        # fmt=".",
        # label="{} messages".format(numm),
        # ms=3
    # )
plt.legend()
plt.title("")
plt.xlabel("")
plt.ylabel("")
plt.savefig("img/latency.pdf")
