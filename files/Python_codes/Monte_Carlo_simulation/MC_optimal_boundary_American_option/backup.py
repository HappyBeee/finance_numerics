import numpy as np

def sample_paths(r, sigma, S_0, T, M, N):
	paths = S_0 * np.ones((N, 1));
	for i in range(M):
		normal_variables = np.random.normal(0, 1, (N, 1));
		ratios = np.exp((r-0.5*sigma*sigma)*T/M + normal_variables * sigma*(T/M)**0.5);
		paths = np.concatenate((paths, paths[:, -1:]*ratios), axis=1);
	return paths;

def optimal_boundary_eval(S, P, K):
	data = np.array([S, P])
	data = np.transpose(data)
	data = sorted(data, key=lambda x: x[0])

	if data[0][0] >= K:
		return P, data[0][0]

	max_tot = -float("inf")
	cutoff = 0
	tot = 0
	for i in range(len(S)-1, -1, -1):
		tot += data[i][1]
		tot -= max(0, K-data[i][0])
		if tot > max_tot:
			max_tot = tot
			cutoff = data[i][0]

	for i in range(len(S)):
		if S[i] >= cutoff:
			continue
		P[i] = K-S[i]

	return P, cutoff

def MC_calculate_ame_put(r, sigma, S_0, K, T, M, N):
	paths = sample_paths(r, sigma, S_0, T, M, N)

	option_prices = np.maximum(0, K-paths[:, -1])
	a = np.exp(-r*T/M)

	optimal_boundary = np.zeros(M+1)
	optimal_boundary[-1] = K

	for i in range(M-1, 0, -1):
		option_prices *= a
		# Execute the put option if stock price is lower than the cutoff, 
		# cutoff is the last stock price point that we keep the related discounted option price.
		option_prices, optimal_boundary[i] = optimal_boundary_eval(paths[:, i], option_prices, K)

	option_prices *= a
	print(optimal_boundary)

	return np.average(option_prices)

if __name__ == "__main__":
	r, sigma, S_0, K, T = 0.1, 0.4, 50, 60, 5/12.
	M, N = 200, 40000

	Ame_put_price = MC_calculate_ame_put(r, sigma, S_0, K, T, M, N)

	print("steps: {}, paths: {} \nAmerican put option price: {:.5f}".format(M, N, Ame_put_price))
