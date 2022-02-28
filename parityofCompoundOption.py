import numpy as np
import matplotlib.pyplot as plt

from BinomialTree import BinomialTree

if __name__ == '__main__':

    S0 = 343.11         # the closing price of Microsoft on December 29, 2021
    K1 = 25             # strike price of derivative Y
    K2 = 350            # strike price of call option of Mocrosoft
    sigma = 0.2209      # volatility
    r = 0.00157         # the risk-free interest rate
    delta = 0.0067      # the continuous dividend yield
    T1 = 1.0            # the time to maturity of derivative Y
    T2 = 1.0            # the time to maturity of call option of Mocrosoft
    n = 500             # step numbers

    S0s = np.arange(100, 601, 20)           # different stock price

    c_results = np.zeros(len(S0s))          # the price of the call into which the options can be exercised at T1   
    p_results = np.zeros(len(S0s))          # the price of the put into which the options can be exercised at T1 
    cc_results = np.zeros(len(S0s))         # the price of the call on the call
    pc_results = np.zeros(len(S0s))         # the price of the put on the call
    cp_results = np.zeros(len(S0s))         # the price of the call on the put
    pp_results = np.zeros(len(S0s))         #  the price of the call on the put

    for i in range(0, len(S0s)):
        Tree = BinomialTree(S0s[i], K1, T1, K2, T2, sigma, r, delta, n, option='euro')
        cc_results[i], __, __, __, c_results[i] = Tree.call_CallOption_tree()
        pc_results[i], __, __, __, __ = Tree.put_CallOption_tree()
        cp_results[i], __, __, __, p_results[i] = Tree.call_PutOption_tree()
        pp_results[i], __, __, __, __ = Tree.put_PutOption_tree()
        
        print('Completed calculatiion of parity formulas with S0 = {:.2f}.'.format(S0s[i]))

    parity1_left = cc_results - pc_results
    parity1_right = c_results - K1 * np.exp(-r * T1)

    parity2_left = cp_results - pp_results
    parity2_right= p_results - K1 * np.exp(-r * T1)

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0s, c_results, 'b:', marker='o', label='c value', alpha=0.9)
    ax.plot(S0s, p_results, 'r:', marker='d', label='p value', alpha=0.9)
    ax.plot(S0s, cc_results, 'g:', marker='^', label='cc value', alpha=0.9)
    ax.plot(S0s, pc_results, 'y:', marker='<', label='pc value', alpha=0.9)
    ax.plot(S0s, cp_results, 'c:', marker='v', label='cp value', alpha=0.9)
    ax.plot(S0s, pp_results, 'm:', marker='>', label='pp value', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Option price', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('option price of compound option.png', dpi=400, bbox_inches='tight')
        

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(S0s, parity1_left, 'b:', marker='^', label='cc-pc', alpha=0.9)
    ax.plot(S0s, parity1_right, 'r:', marker='v', label='c-${K_1e^{-rT_1}}$', alpha=0.9)
    ax.plot(S0s, parity2_left, 'g:', marker='D', label='cp-pp', alpha=0.9)
    ax.plot(S0s, parity2_right, 'm:', marker='o', label='p-${K_1e^{-rT_1}}$', alpha=0.9)
    ax.set_xlabel('Stock price(\$)', fontsize=14)
    ax.set_ylabel('Parity value', fontsize=14)
    ax.set_xlim([100, 600])
    # ax.set_ylim([10, 70])
    ax.tick_params(direction='in', grid_alpha=0.5)
    ax.legend(loc='best')
    plt.grid(linestyle='-.')
    plt.savefig('parityFomula.png', dpi=400, bbox_inches='tight')

    plt.show()