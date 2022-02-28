# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from datetime import date

from BinomialTree import BinomialTree

if __name__ == '__main__':

    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield

    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft

    n = 500             # step number
    option = 'euro'     # 'euro'——European Options, 'amer'——American Options
    # option = 'amer'     # 'euro'——European Options, 'amer'——American Options


# ----------------------------------------------------------------------------------------------------------------------------
    # Greeks delta, gamma and theta for Derivative Y in the European case

    Time = [3.0 / 12, 1.0, 3.0]
    S0_d = np.arange(100, 601, 10)

    price_euro = np.zeros((len(Time), len(S0_d)))
    Delta_euro = np.zeros((len(Time), len(S0_d)))
    gamma_euro = np.zeros((len(Time), len(S0_d)))
    theta_euro = np.zeros((len(Time), len(S0_d)))

    price_amer = np.zeros((len(Time), len(S0_d)))
    Delta_amer = np.zeros((len(Time), len(S0_d)))
    gamma_amer = np.zeros((len(Time), len(S0_d)))
    theta_amer = np.zeros((len(Time), len(S0_d)))

    for i in range(0, len(Time)):
        T1 = T2 = Time[i]
        for j in range(0, len(S0_d)):
            Tree = BinomialTree(S0_d[j], K1, T1, K2, T2, sigma, r, delta, n, option='euro')
            price_euro[i, j], Delta_euro[i, j], gamma_euro[i, j], theta_euro[i, j], __ = Tree.put_CallOption_tree()

            Tree = BinomialTree(S0_d[j], K1, T1, K2, T2, sigma, r, delta, n, option='amer')
            price_amer[i, j], Delta_amer[i, j], gamma_amer[i, j], theta_amer[i, j], __ = Tree.put_CallOption_tree()

            print('Completed calculation of delta, gamma and theta with T_1 = T_2 = {:.2f}, S0 = {:.2f}.'.format(Time[i], S0_d[j]))

    T1 = T2 = 1.0

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, Delta_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, Delta_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, Delta_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Delta', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroDelta.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, Delta_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, Delta_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, Delta_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Delta', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerDelta.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, gamma_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, gamma_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, gamma_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Gamma', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroGamma.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, gamma_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, gamma_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, gamma_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Gamma', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerGamma.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, theta_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, theta_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, theta_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Theta', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroTheta.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_d, theta_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_d, theta_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_d, theta_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Theta', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerTheta.png', dpi=400, bbox_inches='tight')


# ----------------------------------------------------------------------------------------------------------------------------
#     # Greeks Vega for Derivative Y in the European case
    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield

    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft

    n = 500             # step number

    Time = [3.0 / 12, 1.0, 3.0]
    S0_v = np.arange(100, 601, 10)

    euro_result0 = np.zeros((len(Time), len(S0_v)))
    euro_result1 = np.zeros((len(Time), len(S0_v)))

    amer_result0 = np.zeros((len(Time), len(S0_v)))
    amer_result1 = np.zeros((len(Time), len(S0_v)))

    for i in range(0, len(Time)):
        T1 = T2 = Time[i]
        for j in range(0, len(S0_v)):
            Tree_0 = BinomialTree(S0_v[j], K1, T1, K2, T2, sigma, r, delta, n, option='euro')
            euro_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_v[j]), K1, T1, K2, T2, sigma+0.01, r, delta, n, option='euro')
            euro_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            Tree_0 = BinomialTree(S0_v[j], K1, T1, K2, T2, sigma, r, delta, n, option='amer')
            amer_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_v[j]), K1, T1, K2, T2, sigma + 0.01, r, delta, n, option='amer')
            amer_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            print('Completed calculation of Vega with T_1 = T_2 = {:.2f}, S0 = {:.2f}.'.format(Time[i], S0_v[j]))


    vega_euro = euro_result1 - euro_result0
    vega_amer = amer_result1 - amer_result0

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_v, vega_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_v, vega_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_v, vega_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Vega', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroVega.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_v, vega_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_v, vega_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_v, vega_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Vega', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerVega.png', dpi=400, bbox_inches='tight')


# ----------------------------------------------------------------------------------------------------------------------------
#     # Greeks Rho for Derivative Y in the European case
    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield

    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft

    n = 500             # step number

    Time = [3.0 / 12, 1.0, 3.0]
    S0_r = np.arange(100, 601, 20)

    euro_result0 = np.zeros((len(Time), len(S0_r)))
    euro_result1 = np.zeros((len(Time), len(S0_r)))

    amer_result0 = np.zeros((len(Time), len(S0_r)))
    amer_result1 = np.zeros((len(Time), len(S0_r)))

    for i in range(0, len(Time)):
        T1 = T2 = Time[i]
        for j in range(0, len(S0_r)):
            Tree_0 = BinomialTree(S0_r[j], K1, T1, K2, T2, sigma, r, delta, n, option='euro')
            euro_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_r[j]), K1, T1, K2, T2, sigma, r + 0.01, delta, n, option='euro')
            euro_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            Tree_0 = BinomialTree(S0_r[j], K1, T1, K2, T2, sigma, r, delta, n, option='amer')
            amer_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_r[j]), K1, T1, K2, T2, sigma, r + 0.01, delta, n, option='amer')
            amer_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            print('Completed calculation of Rho with T_1 = T_2 = {:.2f}, S0 = {:.2f}.'.format(Time[i], S0_r[j]))

    rho_euro = euro_result1 - euro_result0
    rho_amer = amer_result1 - amer_result0

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_r, rho_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_r, rho_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_r, rho_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Rho', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroRho.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_r, rho_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_r, rho_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_r, rho_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Rho', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerRho.png', dpi=400, bbox_inches='tight')


# ----------------------------------------------------------------------------------------------------------------------------
#     # Greeks Psi for Derivative Y in the European case
    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield

    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft

    n = 500             # step number

    Time = [3.0 / 12, 1.0, 3.0]
    S0_p = np.arange(100, 601, 20)

    euro_result0 = np.zeros((len(Time), len(S0_p)))
    euro_result1 = np.zeros((len(Time), len(S0_p)))

    amer_result0 = np.zeros((len(Time), len(S0_p)))
    amer_result1 = np.zeros((len(Time), len(S0_p)))

    for i in range(0, len(Time)):
        T1 = T2 = Time[i]
        for j in range(0, len(S0_p)):
            Tree_0 = BinomialTree(S0_p[j], K1, T1, K2, T2, sigma, r, delta, n, option='euro')
            euro_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_p[j]), K1, T1, K2, T2, sigma, r, delta + 0.01, n, option='euro')
            euro_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            Tree_0 = BinomialTree(S0_p[j], K1, T1, K2, T2, sigma, r, delta, n, option='amer')
            amer_result0[i, j] = Tree_0.put_CallOption_tree()[0]

            Tree_1 = BinomialTree((S0_p[j]), K1, T1, K2, T2, sigma, r, delta + 0.01, n, option='amer')
            amer_result1[i, j] = Tree_1.put_CallOption_tree()[0]

            print('Completed calculatiion of Psi with T_1 = T_2 = {:.2f}, S0 = {:.2f}.'.format(Time[i], S0_p[j]))


    psi_euro = euro_result1 - euro_result0
    psi_amer = amer_result1 - amer_result0

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_p, psi_euro[0, :], 'b-', marker='o', label='European case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_p, psi_euro[1, :], 'r-', marker='^', label='European case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_p, psi_euro[2, :], 'k-', marker='+', label='European case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Psi', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('euroPsi.png', dpi=400, bbox_inches='tight')

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0_p, psi_amer[0, :], 'b-', marker='o', label='American case with ${T_1=T_2}$=3 months', alpha=0.9)
    ax.plot(S0_p, psi_amer[1, :], 'r-', marker='^', label='American case with ${T_1=T_2}$=1 year', alpha=0.9)
    ax.plot(S0_p, psi_amer[2, :], 'k-', marker='+', label='American case with ${T_1=T_2}$=3 years', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Psi', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('amerPsi.png', dpi=400, bbox_inches='tight')


# ----------------------------------------------------------------------------------------------

    plt.show()
