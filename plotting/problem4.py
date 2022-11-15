# !python
# -*- coding: utf-8 -*

__author__ = 'Erling Ween Eriksen'
__email__ = 'erlinge@nmbu.no'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

MECX4 = pd.read_csv(r'..\MECX_4.csv', names=['M', 'E']).astype({'M': 'float32', 'E': 'float32'})

k = 1.38 * 10 ** (-23)
T = 1 / k
E_ana = -1.962
E2_ana = 15.5558
M_ana = 0.981
M2_ana = 2.9237
C_ana = 4.303 * 10 ** (-26)
X_ana = 0.73162

n_iter = len(MECX4.index)
plot_n = 1000000 + 1
plot_range = np.arange(0, plot_n)
N_spins = 4

E_exp = np.cumsum(MECX4.E.values)/(plot_range + 1)
E2_exp = np.cumsum(MECX4.E.values ** 2) / (plot_range + 1)
M_exp = np.cumsum(MECX4.E.values) / (plot_range + 1)
M2_exp = np.cumsum(MECX4.E.values ** 2) / (plot_range + 1)
heatCapacity = (1 / (N_spins * T * T)) * (E2_exp - E_exp ** 2)

# E_exp = [MECX4.E.iloc[0:ind].sum()/(ind + 1) for ind in MECX4.index]
# M_exp = [MECX4.M.iloc[0:ind].sum()/(ind + 1) for ind in MECX4.index]


plt.plot(E_exp[0:plot_n], label='E Markov Monte Carlo')
plt.plot([E_ana for _ in plot_range], label='E analytical')
plt.legend()
plt.show()
plt.plot(M_exp[0:plot_n], label='M Markov Monte Carlo')
plt.plot([M_ana for _ in plot_range], label='M analytical')
plt.legend()
plt.show()

plt.plot(MECX4.E2.iloc[0:plot_n], label='E2 Markov Monte Carlo')
plt.plot([E2_ana for i in plot_range], label='E2 analytical')

plt.legend()
plt.show()
plt.plot(MECX4.M2.iloc[0:plot_n], label='M2 Markov Monte Carlo')
plt.plot([M2_ana for _ in plot_range], label='M2 analytical')
plt.legend()
plt.show()

plt.plot(MECX4.C.iloc[0:plot_n], label='C Markov Monte Carlo')
plt.plot([C_ana for _ in range(plot_n)], label='C analytical')
plt.legend()
plt.show()
plt.plot(MECX4.X.iloc[0:plot_n], label='X Markov Monte Carlo')
plt.plot([X_ana for _ in range(plot_n)], label='X analytical')
plt.legend()
plt.show()
