# Valuation of European Call Option
# via Monte Carlo Simulation

import numpy as np
from scipy.optimize import fsolve
import scipy.stats as st
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt

from BinomialTree import BinomialTree
from MCS import MonteCarlo
from Geske import Geske

if __name__ == '__main__':

    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield

    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft
    num = 1000000

    def eu_BS(S0_est):
        d_1 = (np.log(S0_est / K2) + (r - delta + 1 / 2 * sigma**2) * T2) / (sigma * np.sqrt(T2))
        d_2 = d_1 - sigma * np.sqrt(T2)
        cvalue = S0_est * np.exp(-delta * T2) * st.norm.cdf(d_1) \
            - K2 * np.exp(-r * T2) * st.norm.cdf(d_2)
        return cvalue - K1

    Time = [1.0]
    S0s = np.arange(100, 601, 20)

    S_1 = fsolve(eu_BS, 300)

    MCS_results = np.zeros((len(Time), len(S0s)))
    BT_results = np.zeros((len(Time), len(S0s)))
    Geske_result = np.zeros((len(Time), len(S0s)))
    for i in range(0, len(Time)):
        T1 = T2 = Time[i]
        for j in range(0, len(S0s)):  
            MCS = MonteCarlo(S0s[j], K1, T1, K2, T2, sigma, r, delta, num)
            __, __, __, MCS_results[i, j] = MCS.euroCallbyMCS()

            Tree = BinomialTree(S0s[j], K1, T1, K2, T2, sigma, r, delta, n=500, option='euro')              
            BT_results[i, j], __, __, __, __ = Tree.call_CallOption_tree()
            ges = Geske(S0s[j], S_1, K1, T1, K2, T2, sigma, r, delta)
            Geske_result[i, j] = ges.geske_formula()
            print('Completed calculatiion of result with S0 = {:.2f}.'.format(S0s[j]))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0s, MCS_results[0, :], 'b--', marker='o', label='Monte Carlo Method', alpha=0.9)
    ax.plot(S0s, BT_results[0, :], 'r--', marker='^', label='Binomial Method', alpha=0.9)
    ax.plot(S0s, Geske_result[0, :], 'g--', marker='d', label='Geske Method', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Price of Derivative Y', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('Price of Derivative Y by MCS.png', dpi=400, bbox_inches='tight')

    plt.show()
