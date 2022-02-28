# Valuation of European Call Option
# via Monte Carlo Simulation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats as st
from scipy.stats import multivariate_normal


class Geske(object):
    def __init__(self, S0, S_1, K1, T1, K2, T2, sigma, r, delta):
        self.S0 = float(S0)
        self.S_1 = float(S_1)
        self.K1 = float(K1)
        self.K2 = float(K2)
        self.T1 = float(T1)
        self.T2 = float(T2)
        self.sigma = float(sigma)
        self.r = float(r)
        self.delta = float(delta)

    def geske_formula(self):
        # calculate price of call on call  with Geske formula
        a_1 = (np.log(self.S0 / self.S_1) + (self.r - self.delta + 1 / 2 * self.sigma**2) * self.T1) / (self.sigma * np.sqrt(self.T1))
        a_2 = a_1 - self.sigma * np.sqrt(self.T1)
        b_1 = (np.log(self.S0 / self.K2) + (self.r - self.delta + 1 / 2 * self.sigma**2) * (self.T1 + self.T2)) / (self.sigma * np.sqrt(self.T1 + self.T2))
        b_2 = b_1 - self.sigma * np.sqrt(self.T1 + self.T2)

        mult = multivariate_normal([0, 0], [[1, np.sqrt(self.T1 / (self.T1 + self.T2))], [np.sqrt(self.T1 / (self.T1 + self.T2)), 1]])

        geskeresult = self.S0 * np.exp(-self.delta * (self.T1 + self.T2)) * mult.cdf((a_1, b_1)) -\
            self.K2 * np.exp(-self.r * (self.T1 + self.T2)) * mult.cdf((a_2, b_2)) -\
            self.K1 * np.exp(-self.r * self.T1) * st.norm.cdf(a_2)

        return geskeresult
