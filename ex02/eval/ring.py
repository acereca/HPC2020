import pandas as pd
import glob
import argparse as ap
import matplotlib.pyplot as plt

p = ap.ArgumentParser()
p.add_argument("slurm_job_id", metavar="J", type=int, nargs=1)
args = p.parse_args()

df = pd.DataFrame()
for f in glob.glob("out/" + str(args.slurm_job_id[0]) + "_*"):
    df = df.append(pd.read_csv(f, sep=";"))

plt.errorbar(
    df['nproc'].unique(),
    [df[df['nproc'] == s]['t_msg/ns'].mean() for s in df['nproc'].unique()],
    fmt="."
)
plt.show()
