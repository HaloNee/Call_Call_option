# Valuation of European Call Option
# via Monte Carlo Simulation

import numpy as np
import matplotlib.pyplot as plt

class MonteCarlo:
    def __init__(self, S0, K1, T1, K2, T2, sigma, r, delta, num):
        self.S0 = float(S0)
        self.K1 = float(K1)
        self.K2 = float(K2)
        self.sigma = float(sigma)
        self.r = float(r)
        self.delta = float(delta)
        self.T1 = float(T1)
        self.T2 = float(T2)
        self.T = self.T1 + self.T2
        self.num = num

    def euroCallbyMCS(self):
        rand = np.random.standard_normal(self.num)                      # generate pseude-random numbers
        ST = self.S0 * np.exp((self.r - self.delta - 0.5 * self.sigma ** 2) * (self.T1 + self.T2) +
                              self.sigma * np.sqrt(self.T1 + self.T2) * rand)      # stock price at T1+T2
        pv1 = np.exp(-self.r * self.T2) * np.maximum(ST - self.K2, 0)   #  Call option discount to time T1
        pv0 = np.exp(-self.r * self.T1) * np.maximum(pv1 - self.K1, 0)  # Derivative Y at time 0
        pv = np.sum(pv0) / self.num                                     # avarage price
        return ST, pv1, pv0, pv
