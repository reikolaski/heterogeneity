#%% import
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes', labelsize=18)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=10)    # legend fontsize

#%% plot x**2
f = lambda x: x**2
x = np.linspace(-5, 5, 100)

fig, ax = plt.subplots()
ax.plot(x, f(x), label='$x^2$')
ax.set_xlabel('$x$')
ax.set_ylabel('$x^2$')
ax.legend(frameon=False)
ax.set_ylim(0, None)
fig.savefig('x-squared-py.pdf', bbox_inches='tight')

# %% minimize x**2
res = minimize(f, x0=5)
print(f'Minimizer: {np.round(res.x[0], 3)}')
print(f'Minimum: {np.round(f(res.x[0]), 3)}')

# %%
