import numpy as np
import MC_OLS_American_option as MC
import matplotlib.pyplot as plt

r, sigma, S_0, K, T = 0.1, 0.4, 50, 50, 5/12.
steps, paths = 100, 1000

data = MC.sample_paths(r, sigma, S_0, T, steps, paths)

option_prices = np.maximum(K-data[:, -1], 0)

for t_n in range(steps-1, int(steps/2), -1):
	option_prices *= np.exp(-r*T/steps)
	option_prices = MC.linear_fitting(data[:, t_n], option_prices)
	option_prices = np.maximum(option_prices, K-data[:, t_n])

X = data[:, int(steps/2)]
Y = option_prices*np.exp(-r*T/steps)
Y_fitted = MC.linear_fitting(X, Y)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.scatter(X, Y, c="blue", s=2)
ax.scatter(X, Y_fitted, c="red", s=2)
ax.set_xlabel("stock prices")
ax.set_ylabel("option prices")
ax.set_title("option prices fitting performance")


fig.savefig("linear_fitting.png")
fig.show()