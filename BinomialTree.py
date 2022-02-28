# Valuation of European Call Option
# via Binomial Tree method
import numpy as np
import scipy.stats as st


class BinomialTree:
    def __init__(self, S0, K1, T1, K2, T2, sigma, r, delta, n, option):
        self.S0 = float(S0)
        self.K1 = float(K1)
        self.K2 = float(K2)
        self.sigma = float(sigma)
        self.r = float(r)
        self.delta = float(delta)
        self.T1 = float(T1)
        self.T2 = float(T2)
        self.n = n
        self.h = (self.T1 + self.T2) / n
        self.u = np.exp((self.r - self.delta) * self.h + self.sigma * np.sqrt(self.h))
        self.d = np.exp((self.r - self.delta) * self.h - self.sigma * np.sqrt(self.h))
        self.p = (np.exp((self.r - self.delta) * self.h) - self.d) / (self.u - self.d)
        self.option = option  # 'euro'——European Options, 'amer'——American Options

    def Pro(self, n, j, p):
        # calculate probability for binomial distribution with parameters n, j, p
        pro = 1.0
        num1 = j
        num2 = n - j
        for i in range(n, 0, -1):
            if num1 > num2:
                pro = pro * i / num1 * p
                num1 -= 1
            else:
                pro = pro * i / num2 * (1 - p)
                num2 -= 1
        return pro

    def euroCallOptionTree(self, S_ex):
        # calculate call option price with Binomial Tree method
        c_n = int(self.T2 / self.h)
        s_prices = np.zeros(c_n + 1)        # to store the stock price at time (i)
        c_values = np.zeros(c_n + 1)      # to store the price of call option at time (i)
        result = 0

        k = int(np.ceil((np.log(self.K2 / S_ex) - np.log(self.d) * c_n) / np.log(self.u / self.d)))

        for i in range(k, c_n + 1, 1):
            s_prices[i] = S_ex * self.u**i * self.d ** (c_n - i)
            c_values[i] = np.maximum(s_prices[i] - self.K2, 0)
            result += c_values[i] * self.Pro(c_n, i, self.p)
        return result

    def euroPutOptionTree(self, S_ex):
        # calculate put option price with Binomial Tree method
        c_n = int(self.T2 / self.h)
        s_prices = np.zeros(c_n + 1)        # to store the stock price at time (i)
        c_values = np.zeros(c_n + 1)      # to store the price of call option at time (i)
        result = 0

        k = int(np.floor((np.log(self.K2 / S_ex) - np.log(self.d) * c_n) / np.log(self.u / self.d)))

        for i in range(0, k, 1):
            s_prices[i] = S_ex * self.u**i * self.d ** (c_n - i)
            c_values[i] = np.maximum(self.K2 - s_prices[i], 0)
            result += c_values[i] * self.Pro(c_n, i, self.p)
        return result

    def call_BS_formula(self, S_ex):
        # calculate call option price with B-S formula
        d_1 = (np.log(S_ex / self.K2) + (self.r - self.delta + 1 / 2 * self.sigma**2) * self.T2) / (self.sigma * np.sqrt(self.T2))
        d_2 = d_1 - self.sigma * np.sqrt(self.T2)
        callOptionPrice = S_ex * np.exp(-self.delta * self.T2) * st.norm.cdf(d_1) \
            - self.K2 * np.exp(-self.r * self.T2) * st.norm.cdf(d_2)

        return callOptionPrice

    def put_BS_formula(self, S_ex):
        # calculate put option price with B-S formula
        d_1 = (np.log(S_ex / self.K2) + (self.r - self.delta + 1 / 2 * self.sigma**2) * self.T2) / (self.sigma * np.sqrt(self.T2))
        d_2 = d_1 - self.sigma * np.sqrt(self.T2)
        putOptionPrice = self.K2 * np.exp(-self.r * self.T2) * st.norm.cdf(-d_2) \
            - S_ex * np.exp(-self.delta * self.T2) * st.norm.cdf(-d_1)

        return putOptionPrice

    def call_CallOption_tree(self):
        prices = np.zeros(self.n + 1)        # to store the stock price at time (i)
        c_values = np.zeros(self.n + 1)      # to store the price of call option at time (i)
        cc_values = np.zeros((self.n + 1, self.n + 1))     # to store the value of Derivative Y at time (i)
        for i in range(0, self.n + 1):
            prices[i] = self.S0 * self.u**i * self.d ** (self.n - i)
            c_values[i] = np.maximum(prices[i] - self.K2, 0)
        for j in range(self.n - 1, -1, -1):
            for i in range(0, j + 1):
                prices[i] = prices[i + 1] / self.u

                if int(j) > int(self.n * self.T1 / (self.T1 + self.T2)):          # time between (T1, T2)
                    c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])

                elif int(j) == int(self.n * self.T1 / (self.T1 + self.T2)):       # time at T1
                    if self.option == 'euro':
                        c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])
                    elif self.option == 'amer':
                        c_values[i] = self.call_BS_formula(prices[i])
                        # c_values[i] = self.euroCallOptionTree(prices[i])  # another methed using binomial tree

                    cc_values[j, i] = np.maximum(c_values[i] - self.K1, 0)

                else:                                                   # time between(0, T1)
                    if self.option == 'euro':
                        c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])
                        cc_values[j, i] = np.exp(-self.r * self.h) * (self.p * cc_values[j + 1, i + 1] + (1 - self.p) * cc_values[j + 1, i])
                    elif self.option == 'amer':
                        c_values[i] = self.call_BS_formula(prices[i])
                        # c_values[i] = self.euroCallOptionTree(prices[i])  # another methed using binomial tree
                        cc_values[j, i] = np.maximum(c_values[i] - self.K1, np.exp(-self.r * self.h) *
                                                     (self.p * cc_values[j + 1, i + 1] + (1 - self.p) * cc_values[j + 1, i]))

        Delta_0 = np.exp(-self.r * self.h) * (cc_values[1, 1] - cc_values[1, 0]) / (self.u - self.d) / self.S0

        Delta_u = np.exp(-self.r * self.h) * (cc_values[2, 2] - cc_values[2, 1]) / (self.u * self.u - self.u * self.d) / self.S0
        Delta_d = np.exp(-self.r * self.h) * (cc_values[2, 1] - cc_values[2, 0]) / (self.u * self.d - self.d * self.d) / self.S0
        gamma_0 = (Delta_u - Delta_d) / (self.u - self.d) / self.S0

        epsilon = self.u * self.d * self.S0 - self.S0
        theta_0 = (cc_values[2, 1] - epsilon * Delta_0 - 1 / 2 * epsilon**2 * gamma_0 - cc_values[0, 0]) / (2 * self.h)
        theta_0 = theta_0 / 365

        return cc_values[0, 0], Delta_0, gamma_0, theta_0, c_values[0]

    def put_CallOption_tree(self):
        prices = np.zeros(self.n + 1)        # to store the stock price at time (i)
        c_values = np.zeros(self.n + 1)      # to store the price of call option at time (i)
        pc_values = np.zeros((self.n + 1, self.n + 1))     # to store the value of Derivative Y at time (i)
        for i in range(0, self.n + 1):
            prices[i] = self.S0 * self.u**i * self.d ** (self.n - i)
            c_values[i] = np.maximum(prices[i] - self.K2, 0)
            # print(prices[i], c_values[i])
        for j in range(self.n - 1, -1, -1):
            for i in range(0, j + 1):
                prices[i] = prices[i + 1] / self.u

                if int(j) > int(self.n * self.T1 / (self.T1 + self.T2)):          # time between (T1, T2)
                    c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])

                elif int(j) == int(self.n * self.T1 / (self.T1 + self.T2)):       # time at T1
                    if self.option == 'euro':
                        c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])
                    elif self.option == 'amer':
                        c_values[i] = self.call_BS_formula(prices[i])
                        # c_values[i] = self.euroCallOptionTree(prices[i])  # another methed using binomial tree

                    pc_values[j, i] = np.maximum(self.K1 - c_values[i], 0)

                else:                                                   # time between(0, T1)
                    if self.option == 'euro':
                        c_values[i] = np.exp(-self.r * self.h) * (self.p * c_values[i + 1] + (1 - self.p) * c_values[i])
                        pc_values[j, i] = np.exp(-self.r * self.h) *\
                            (self.p * pc_values[j + 1, i + 1] + (1 - self.p) * pc_values[j + 1, i])
                    elif self.option == 'amer':
                        c_values[i] = self.call_BS_formula(prices[i])
                        # c_values[i] = self.euroCallOptionTree(prices[i])  # another methed using binomial tree
                        pc_values[j, i] = np.maximum(self.K1 - c_values[i], np.exp(-self.r * self.h) *
                                                     (self.p * pc_values[j + 1, i + 1] + (1 - self.p) * pc_values[j + 1, i]))


        Delta_0 = np.exp(-self.r * self.h) * (pc_values[1, 1] - pc_values[1, 0]) / (self.u - self.d) / self.S0

        Delta_u = np.exp(-self.r * self.h) * (pc_values[2, 2] - pc_values[2, 1]) / (self.u * self.u - self.u * self.d) / self.S0
        Delta_d = np.exp(-self.r * self.h) * (pc_values[2, 1] - pc_values[2, 0]) / (self.u * self.d - self.d * self.d) / self.S0
        gamma_0 = (Delta_u - Delta_d) / (self.u - self.d) / self.S0

        epsilon = self.u * self.d * self.S0 - self.S0
        theta_0 = (pc_values[2, 1] - epsilon * Delta_0 - 1 / 2 * epsilon**2 * gamma_0 - pc_values[0, 0]) / (2 * self.h)
        theta_0 = theta_0 / 365

        return pc_values[0, 0], Delta_0, gamma_0, theta_0, c_values[0]

    def call_PutOption_tree(self):
        prices = np.zeros(self.n + 1)        # to store the stock price at time (i)
        p_values = np.zeros(self.n + 1)      # to store the price of call option at time (i)
        cp_values = np.zeros((self.n + 1, self.n + 1))     # to store the value of Derivative Y at time (i)
        for i in range(0, self.n + 1):
            prices[i] = self.S0 * self.u**i * self.d ** (self.n - i)
            p_values[i] = np.maximum(self.K2 - prices[i], 0)
            # print(prices[i], p_values[i])
        for j in range(self.n - 1, -1, -1):
            for i in range(0, j + 1):
                prices[i] = prices[i + 1] / self.u

                if int(j) > int(self.n * self.T1 / (self.T1 + self.T2)):          # time between (T1, T2)
                    p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])

                elif int(j) == int(self.n * self.T1 / (self.T1 + self.T2)):       # time at T1
                    if self.option == 'euro':
                        p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])
                    elif self.option == 'amer':
                        p_values[i] = self.put_BS_formula(prices[i])
                        # p_values[i] = self.euroPutOptionTree(prices[i])  # another methed using binomial tree

                    cp_values[j, i] = np.maximum(p_values[i] - self.K1, 0)

                else:                                                   # time between(0, T1)
                    if self.option == 'euro':
                        p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])
                        cp_values[j, i] = np.exp(-self.r * self.h) * (self.p * cp_values[j + 1, i + 1] + (1 - self.p) * cp_values[j + 1, i])
                    elif self.option == 'amer':
                        p_values[i] = self.put_BS_formula(prices[i])
                        # p_values[i] = self.euroPutOptionTree(prices[i])  # another methed using binomial tree
                        cp_values[j, i] = np.maximum(p_values[i] - self.K1, np.exp(-self.r * self.h) *
                                                     (self.p * cp_values[j + 1, i + 1] + (1 - self.p) * cp_values[j + 1, i]))


        Delta_0 = np.exp(-self.r * self.h) * (cp_values[1, 1] - cp_values[1, 0]) / (self.u - self.d) / self.S0

        Delta_u = np.exp(-self.r * self.h) * (cp_values[2, 2] - cp_values[2, 1]) / (self.u * self.u - self.u * self.d) / self.S0
        Delta_d = np.exp(-self.r * self.h) * (cp_values[2, 1] - cp_values[2, 0]) / (self.u * self.d - self.d * self.d) / self.S0
        gamma_0 = (Delta_u - Delta_d) / (self.u - self.d) / self.S0

        epsilon = self.u * self.d * self.S0 - self.S0
        theta_0 = (cp_values[2, 1] - epsilon * Delta_0 - 1 / 2 * epsilon**2 * gamma_0 - cp_values[0, 0]) / (2 * self.h)
        theta_0 = theta_0 / 365

        return cp_values[0, 0], Delta_0, gamma_0, theta_0, p_values[0]

    def put_PutOption_tree(self):
        prices = np.zeros(self.n + 1)        # to store the stock price at time (i)
        p_values = np.zeros(self.n + 1)      # to store the price of call option at time (i)
        pp_values = np.zeros((self.n + 1, self.n + 1))     # to store the value of Derivative Y at time (i)
        for i in range(0, self.n + 1):
            prices[i] = self.S0 * self.u**i * self.d ** (self.n - i)
            p_values[i] = np.maximum(self.K2 - prices[i], 0)
            # print(prices[i], p_values[i])
        for j in range(self.n - 1, -1, -1):
            for i in range(0, j + 1):
                prices[i] = prices[i + 1] / self.u

                if int(j) > int(self.n * self.T1 / (self.T1 + self.T2)):          # time between (T1, T2)
                    p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])

                elif int(j) == int(self.n * self.T1 / (self.T1 + self.T2)):       # time at T1
                    if self.option == 'euro':
                        p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])
                    elif self.option == 'amer':
                        p_values[i] = self.put_BS_formula(prices[i])

                        # p_values[i] = self.euroPutOptionTree(prices[i])  # another methed using binomial tree

                    pp_values[j, i] = np.maximum(self.K1 - p_values[i], 0)

                else:                                                   # time between(0, T1)
                    if self.option == 'euro':
                        p_values[i] = np.exp(-self.r * self.h) * (self.p * p_values[i + 1] + (1 - self.p) * p_values[i])
                        pp_values[j, i] = np.exp(-self.r * self.h) * (self.p * pp_values[j + 1, i + 1] + (1 - self.p) * pp_values[j + 1, i])
                    elif self.option == 'amer':
                        p_values[i] = self.put_BS_formula(prices[i])
                        # p_values[i] = self.euroPutOptionTree(prices[i])  # another methed using binomial tree
                        pp_values[j, i] = np.maximum(self.K1 - p_values[i], np.exp(-self.r * self.h) *
                                                     (self.p * pp_values[j + 1, i + 1] + (1 - self.p) * pp_values[j + 1, i]))

        Delta_0 = np.exp(-self.r * self.h) * (pp_values[1, 1] - pp_values[1, 0]) / (self.u - self.d) / self.S0

        Delta_u = np.exp(-self.r * self.h) * (pp_values[2, 2] - pp_values[2, 1]) / (self.u * self.u - self.u * self.d) / self.S0
        Delta_d = np.exp(-self.r * self.h) * (pp_values[2, 1] - pp_values[2, 0]) / (self.u * self.d - self.d * self.d) / self.S0
        gamma_0 = (Delta_u - Delta_d) / (self.u - self.d) / self.S0

        epsilon = self.u * self.d * self.S0 - self.S0
        theta_0 = (pp_values[2, 1] - epsilon * Delta_0 - 1 / 2 * epsilon**2 * gamma_0 - pp_values[0, 0]) / (2 * self.h)
        theta_0 = theta_0 / 365

        return pp_values[0, 0], Delta_0, gamma_0, theta_0, p_values[0]
