// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

def sample_paths(r, sigma, S_0, T, M, N):
	""" Generate N stock price change paths, each has length M+1.
	"""
	paths = list()
	for i in range(N):
		path = [S_0]
		for j in range(M):
			new_price = path[-1]*np.exp((r-0.5*sigma*sigma)*T/M + np.random.normal()*sigma*np.sqrt(T/M))
			path.append(new_price)
		paths.append(path)
	return paths

def optimal_boundary_eval(S, P, K):
	""" Input all stock prices S on every paths at some time, 
		, all reference option prices P, and exercise price K.
		Find the optimal execution stock price as a cutoff.
		Return the option prices after applied the cutoff, and the cutoff value.
		The option prices has the same order with input stock prices.
	""" 
	length = len(S)
	data_pairs = [[S[i], P[i]] for i in range(length)]
	data_pairs = sorted(data_pairs, key=lambda x: x[0])
	new_P = list()

	total = 0
	total_max = 0
	cutoff = 0
	for i in range(length-1, -1, -1):
		total += data_pairs[i][1]
		total -= max(0, K-data_pairs[i][0])
		if total > total_max:
			total_max = total
			cutoff = data_pairs[i][0]
	
	for i in range(length):
		if S[i] < cutoff:
			new_P.append(max(K-S[i], 0))
		else:
			new_P.append(P[i])

	return (new_P, cutoff)

def MC_optimal_boundary_Ame_put(r, sigma, S_0, K, T, M, N):
	""" Main function for calculating American put option with Monte Carlo simulaiton.
	""" 
	put_price = 0
	paths = sample_paths(r, sigma, S_0, T, M, N)
	
	# Stock prices and option prices at exercise time.
	stock_prices = [paths[i][-1] for i in range(N)]
	put_prices = [max(K-stock_prices[i], 0) for i in range(N)]
	
	# Initialize the cutoffs, to store the optimal exercise stock prices, for later on debug.
	cutoffs = [0]*(M+1)
	cutoffs[-1] = K

	# Backward calcuation.
	for i in range(M-1, -1, -1):
		for j in range(N):
			stock_prices[j] = paths[j][i]
			put_prices[j] = put_prices[j] * np.exp(-r*T/M)
		(put_prices, cutoffs[i]) = optimal_boundary_eval(stock_prices, put_prices, K)
	
	# Get the average value of the put prices at the beginning of each paths.
	for i in range(N):
		put_price += put_prices[i]/N

	return put_price

if __name__ == "__main__":
	r, sigma, S_0, K, T = 0.1, 0.4, 50, 60, 5.0/12
	M, N = 100, 10000

	Ame_put_price = MC_optimal_boundary_Ame_put(r, sigma, S_0, K, T, M, N)
	print("MC American put price: {0:.5f}".format(Ame_put_price))
