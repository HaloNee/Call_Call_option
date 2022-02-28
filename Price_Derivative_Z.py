# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from BinomialTree import BinomialTree

if __name__ == '__main__':

    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of Derivative Z
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209      # volatility
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield
    T1 = 1.0            # the time to maturity of Derivative Z
    T2 = 1.0            # the time to maturity of call option of Mocrosoft

    n = 500             # step number
    # option = 'euro'     # 'euro'——European Options, 'amer'——American Options
    option = 'amer'     # 'euro'——European Options, 'amer'——American Options

# ---------------------------------------------------------------------------------------------
    # fix(K_1, T_1, K_2, T_2) = (25, 1, 350, 1) and change n from 2 to 200
    steps = np.arange(2, 201, 2)

    results_euro = np.zeros(len(steps))
    results_amer = np.zeros(len(steps))
    for i in np.arange(0, len(steps), 1):
        step = steps[i]
        Tree = BinomialTree(S0, K1, T1, K2, T2, sigma, r, delta, step, option='euro')
        results_euro[i] = Tree.put_CallOption_tree()[0]
        Tree = BinomialTree(S0, K1, T1, K2, T2, sigma, r, delta, step, option='amer')
        results_amer[i] = Tree.put_CallOption_tree()[0]
        print('Completed calculation with step = {:.2f}.'.format(step))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(steps, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    # ax.plot(steps, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('Step number', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([0, 200])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroPrice_on_different_step.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    # ax.plot(steps, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    ax.plot(steps, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('Step number', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([0, 200])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    plt.grid(linestyle='-.')
    plt.savefig('amerPrice_on_different_step.png', dpi=400, bbox_inches='tight')

# --------------------------------------------------------------------------------------------
    # fix(K_1, T_1, K_2, n) = (25, 1, 350, 500) and change T_2 from 0.1 to 3.0
    T2s = np.arange(0.1, 3.1, 0.1)

    results_euro = np.zeros(len(T2s))
    results_amer = np.zeros(len(T2s))
    for i in np.arange(0, len(T2s), 1):
        Tree = BinomialTree(S0, K1, T1, K2, T2s[i], sigma, r, delta, n, option='euro')
        results_euro[i] = Tree.put_CallOption_tree()[0]
        Tree = BinomialTree(S0, K1, T1, K2, T2s[i], sigma, r, delta, n, option='amer')
        results_amer[i] = Tree.put_CallOption_tree()[0]
        print('Completed calculation with T_2 = {:.2f}.'.format(T2s[i]))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(T2s, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    ax.plot(T2s, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('${T_2}$/year', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([0.0, 3.0])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')

    plt.savefig('Price_on_different_T2.png', dpi=400, bbox_inches='tight')

# --------------------------------------------------------------------------------------------
    # fix(K_1, K_2, T_2, n) = (25, 350, 1, 500) and change T_1 from 0.1 to 3.0
    T1s = np.arange(0.1, 3.1, 0.1)
    results_euro = np.zeros(len(T1s))
    results_amer = np.zeros(len(T1s))
    for i in np.arange(0, len(T1s), 1):
        Tree = BinomialTree(S0, K1, T1s[i], K2, T2, sigma, r, delta, n, option='euro')
        results_euro[i] = Tree.put_CallOption_tree()[0]
        Tree = BinomialTree(S0, K1, T1s[i], K2, T2, sigma, r, delta, n, option='amer')
        results_amer[i] = Tree.put_CallOption_tree()[0]
        print('Completed calculation  with T_1 = {:.2f}.'.format(T1s[i]))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(T1s, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    ax.plot(T1s, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('${T_1}$/year', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([0.0, 3.0])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')

    plt.savefig('Price_on_different_T1.png', dpi=400, bbox_inches='tight')

# --------------------------------------------------------------------------------------------
    # fix(K_1, T_1, T_2, n) = (25, 1, 1, 500) and change K_2 from 200 to 500
    K2s = np.arange(200, 501, 10)
    results_euro = np.zeros(len(K2s))
    results_amer = np.zeros(len(K2s))
    for i in np.arange(0, len(K2s), 1):
        Tree = BinomialTree(S0, K1, T1, K2s[i], T2, sigma, r, delta, n, option='euro')
        results_euro[i] = Tree.put_CallOption_tree()[0]

        Tree = BinomialTree(S0, K1, T1, K2s[i], T2, sigma, r, delta, n, option='amer')
        results_amer[i] = Tree.put_CallOption_tree()[0]
        print('Completed calculation with K2 = {:.2f}.'.format(K2))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(K2s, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    ax.plot(K2s, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('${K_2}$(\$)', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([200, 500])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')

    plt.savefig('Price_on_different_K2.png', dpi=400, bbox_inches='tight')

# --------------------------------------------------------------------------------------------
    # fix(T_1, K_2, T_2, n) = (1, 350, 1, 500) and change K_1 from 10 to 50
    K1s = np.arange(10, 51, 2)
    results_euro = np.zeros(len(K1s))
    results_amer = np.zeros(len(K1s))
    for i in np.arange(0, len(K1s), 1):
        Tree = BinomialTree(S0, K1s[i], T1, K2, T2, sigma, r, delta, n, option='euro')
        results_euro[i] = Tree.put_CallOption_tree()[0]

        Tree = BinomialTree(S0, K1s[i], T1, K2, T2, sigma, r, delta, n, option='amer')
        results_amer[i] = Tree.put_CallOption_tree()[0]
        print('Completed calculation with K1 = {:.2f}.'.format(K1))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(K1s, results_euro, 'b-.', marker='o', label='European', alpha=0.9)
    ax.plot(K1s, results_amer, 'r-.', marker='^', label='American', alpha=0.9)
    ax.set_xlabel('${K_1}$(\$)', fontsize=14)
    ax.set_ylabel('Price of Derivative Z', fontsize=14)
    ax.set_xlim([10, 50])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')

    plt.savefig('Price_on_different_K1.png', dpi=400, bbox_inches='tight')

# --------------------------------------------------------------------------------------------

    plt.show()
