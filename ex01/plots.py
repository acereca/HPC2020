import matplotlib as mpl
from matplotlib import rc as rc
from matplotlib import pyplot as plt

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.style.use('bmh')

import numpy as np

### 1.2
def moores_law(N0, t):
    return N0 * 2**(1/18 * t)

x = np.linspace(0, 26, 1000)
y = moores_law(415530e12, x)

plt.figure()
plt.plot(x, y, label=r"$P_{\textrm{compute}}(t) = N_0 2^{\frac{1}{18} t}$")
plt.hlines(1e18, 0, 26, linestyle='--', color='k')
plt.title("Time until exa flop (from June 2020)")
plt.xlabel("time [months]")
plt.ylabel("Performance [flop/s]")
plt.legend()
plt.savefig("./plots/Moore.pdf")

### 1.3
import datetime as dt
from matplotlib import dates as mdt
Rmax_2011 = 10510e12
Rmax_2007 = 478.2e12
x_data = [dt.date(2007,11,1), dt.date(2011,11,1)]
y_data = [Rmax_2007, Rmax_2011]

# get growth rate by linear fit to log values (ignore statistical nonsense due
# to non-weightetness for small values... just an approximation)
params = np.polyfit([2007,2011], [np.log(Rmax_2007), np.log(Rmax_2011)], 1)
poly = np.poly1d(params)

yfit = lambda x: np.exp(poly(x)) # transform values back from log domain
x = [dt.date(2007,11,1),dt.date(2020,11,1)]
y = [yfit(2007), yfit(2020)]

plt.figure()
plt.plot(x, y, label="Linear extrapolation")
plt.plot(x_data, y_data, marker='x', linestyle='none', color='k', label="Rmax [flop/s]")
plt.hlines(1e18,x[0], x[1], linestyle='--', color='k')
plt.yscale('log')
plt.title("Extrapolating from past performance to exa scale")
plt.xlabel("Year")
plt.ylabel("Performance [flop/s]")
plt.legend()
plt.savefig("./plots/GrowthRate.pdf")
