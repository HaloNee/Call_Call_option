import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from datetime import date
from scipy import stats
from scipy.optimize import fsolve


class call_option(object):

	def __init__(self, S0, K, sigma, r, t, M, delta):
		self.S0 = float(S0)
		self.K = K
		self.sigma = sigma
		self.r = r
		self.t = t
		self.M = M
		self.delta = delta


	def update_ttm(self):
		if self.t > self.M:
			raise ValueError("Pricing date later than maturity.")
		self.T = (self.M - self.t).days / 365.

	def d1(self):
		d1 = ((np.log(self.S0 / self.K)
			   + (self.r - self.delta + 0.5 * self.sigma ** 2) * self.T)
			  / (self.sigma * np.sqrt(self.T)))
		return d1

	def d2(self):
		d2 = ((np.log(self.S0 / self.K)
			   + (self.r - self.delta - 0.5 * self.sigma ** 2) * self.T)
			  / (self.sigma * np.sqrt(self.T)))
		return d2


	def value(self):
		''' return present value of call option'''
		self.update_ttm()
		d1 = self.d1()
		d2 = self.d2()
		value = (self.S0 * np.exp(-self.delta * self.T)* stats.norm.cdf(d1, 0.0, 1.0)
				 - self.K * np.exp(-self.r * self.T) * stats.norm.cdf(d2, 0.0, 1.0))
		return value

	def imp_vol(self, C0, sigma_est=0.2):
		''' Return implied volatility given option price. 
		C0 : call option price
		sigma_est : initial estimated value if sigma
		'''
		option = call_option(self.S0, self.K, sigma_est, self.r, self.t,
		 self.M, self.delta)

		option.update_ttm()

		def difference(sigma):
			option.sigma = sigma
			return option.value() - C0
		ivol = fsolve(difference, sigma_est)[0]
		return ivol


if __name__ == '__main__':
	S0 = 343.11		# current share price, at 2021,11,19
	K = 345			
	sigma = 0.2
	r = 0.00157		 
	t = date(2021, 11, 19)
	M = date(2023, 1, 20)
	delta = 0.0067
	sigma_est = 0.2
	
	
	data = pd.read_csv('option priceM=20230120.csv',index_col=0)
	# print(data)
	data['Imp_Vol'] = 0.0
	for row in data.index:
		K = data['Strike'][row]
		C0 = data['Calls'][row]
		call = call_option(S0, K, sigma, r, t, M, delta)
		data.loc[row, 'Imp_Vol'] = call.imp_vol(C0, sigma_est)
	
	fig = plt.figure(figsize=(10, 5))
	ax= fig.add_subplot(1, 1, 1)
	ax.plot(data['Strike'], data['Imp_Vol'], 'k.', alpha=0.9)
	ax.set_xlabel('Strike(\$)',fontsize=14);		# , fontweight='bold'
	ax.set_ylabel('Implied volatilities', fontsize=14)	# , fontweight='bold'
	ax.set_xlim([300, 500])
	ax.set_ylim([0, 0.5])
	ax.tick_params(direction='in', grid_alpha=0.5)
	plt.grid(linestyle='-.')
	
	plt.savefig('Imp_Vol.png', dpi=400, bbox_inches='tight')
	plt.show()

