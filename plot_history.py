from datetime import date
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('wednesdays price.csv',index_col=0, parse_dates=True)

fig = plt.figure(figsize=(9, 6))
ax= fig.add_subplot(1, 1, 1)
data['Close'].plot(ax=ax, alpha=0.9)
ax.set_xlabel('Date',fontsize=14);		# , fontweight='bold'
ax.set_ylabel('Stock price', fontsize=14)	# , fontweight='bold'
ax.set_xlim(['1/1/2011', '1/1/2021'])
ax.set_ylim([0, 250])
ax.tick_params(direction='in', grid_alpha=0.5)
plt.grid(linestyle='-.')

plt.savefig('stock_price.png', dpi=400, bbox_inches='tight')
plt.show()

