import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("ggplot")

dim = "Grid size"
nodes = "# Ranks"
it = "Iterations"
t = "Time to complete / us"
tpit = "Time per iteration / us"
fpit = "Flops per iteration"

data = pd.read_csv(
    "../code/out/grid.csv", sep=";", comment="#", names=[dim, nodes, it, t]
)

# data[t] /= 1e9
data[dim] = data[dim].astype(str) + "x" + data[dim].astype(str)
data[tpit] = data[t] / data[it]
# data[fpit] = 7 * (data['dim'] - 1) * (data['dim'] - 1)
# data['Flops total'] = data['iterations'] * data[fpit]
# data['GFLOP64/s'] = data[fpit] / data[tpit] / 1e9
# data['Grid Size'] = data['dim'].astype(str) + "x" + data['dim'].astype(str)

# data_out = data[['Grid Size', tpit, 'Flops total', 'GFLOP64/s']]

data = data.set_index(dim)


# Times
data_t = pd.DataFrame()

for n in data[nodes].unique():
    data_t['NP = {:02d}'.format(n)] = data[data[nodes] == n][tpit]

data_t = data_t.sort_index(axis=1)

data_t.to_latex(
    "data/heat_t.tex",
    float_format="%d",
    # index=False,
    # escape=False,
    encoding="utf-8",
)

# Speedup
data_s = pd.DataFrame()

for n in data[nodes].unique():
    data_s['NP = {:02d}'.format(n)] = data[data[nodes] == 1][tpit] / data[data[nodes] == n][tpit]

data_s = data_s.sort_index(axis=1)

data_s.to_latex(
    "data/heat_s.tex",
    float_format="%.4f",
    # index=False,
    # escape=False,
    encoding="utf-8",
)

# Efficiency
data_e = pd.DataFrame()

for n in data[nodes].unique():
    data_e['NP = {:02d}'.format(n)] = data[data[nodes] == 1][tpit] / data[data[nodes] == n][tpit] / n

data_e = data_e.sort_index(axis=1)

data_e.to_latex(
    "data/heat_e.tex",
    float_format="%.4f",
    # index=False,
    # escape=False,
    encoding="utf-8",
)

print(data_t)
print(data_s)
print(data_e)
