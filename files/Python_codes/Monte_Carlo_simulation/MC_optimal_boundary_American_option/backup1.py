import numpy as np

def sample_paths(r, sigma, S_0, T, M, N):
	""" 蒙卡抽样股价变化路径。
	"""
	paths = list()
	for i in range(N):
		path = [S_0]
		for j in range(M):
			new_price = path[-1] * np.exp((r-0.5*sigma*sigma)*T/M+sigma*np.sqrt(T/M)*np.random.normal())
			path.append(new_price)
		paths.append(path)
	return paths

def optimal_boundary_eval(S, P, K):
	""" 输入某时刻所有路径上的股票价格和参考期权价格，确定最佳执行边界股价，
		返回考虑执行边界后的期权价格和执行边界。
	"""
	S = list(S)
	P = list(P)
	length = len(S)
	# 使用最佳执行边界处理后的期权价格。
	P_filtered = [None for i in range(length)]

	# 组成股价和对应参考期权价格对，并根据股价由小到大排序。
	data_pairs = [[S[i], P[i]] for i in range(length)]
	data_pairs = sorted(data_pairs, key=lambda x:x[0])

	tot = 0
	max_tot = -float("inf")
	cutoff = 0
	for i in range(length-1, -1, -1):
		tot += data_pairs[i][1]
		tot -= max(K-data_pairs[i][0], 0)
		if tot > max_tot:
			max_tot = tot
			cutoff = data_pairs[i][0]
	for i in range(length):
		if S[i] < cutoff:
			P_filtered[i] = max(K-S[i], 0)
		else:
			P_filtered[i] = P[i]
	return (P_filtered, cutoff)

def MC_optimal_boundary_Ame_put(r, sigma, S_0, K, T, M, N):
	result = 0
	paths = sample_paths(r, sigma, S_0, T, M, N)
	# cutoffs 为执行边界股价，留作debug用。
	cutoffs = [0 for i in range(M+1)]
	cutoffs[-1] = K
	# 执行时刻股价和期权价格。
	stock_prices = [paths[i][-1] for i in range(N)]
	put_prices = [max(K-paths[i][-1],0) for i in range(N)]

	# 由执行时刻往回一步一步递推。
	for i in range(M-1, 0, -1):
		# 确定该时刻股票价格和贴现回来的期权价格。
		for j in range(N):
			stock_prices[j] = paths[j][i]
			put_prices[j] *= np.exp(-r*T/M)
		# 使用最佳执行边界函数重新确定每个股价对应的期权价格。
		put_prices, cutoffs[i] = optimal_boundary_eval(stock_prices, put_prices, K)
	
	# 最后再贴现一步回初始时刻，并计算出平均期权价格。
	for i in range(N):
		put_prices[i] *= np.exp(-r*T/M)
		result += put_prices[i]/N

	print(cutoffs)
	return result

if __name__ == "__main__":
	r, sigma, S_0, K, T = 0.1, 0.4, 50, 60, 5.0/12
	M, N = 100, 10000

	put_price = MC_optimal_boundary_Ame_put(r, sigma, S_0, K, T, M, N)
	print("American put price:{0:.5f}".format(put_price))
