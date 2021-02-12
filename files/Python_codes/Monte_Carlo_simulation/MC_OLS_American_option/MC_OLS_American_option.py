// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

def sample_paths(r, sigma, S_0, T, steps, paths):
	data = S_0 * np.ones((paths, 1))

	for i in range(steps):
		normal_variables = np.random.normal(0, 1, (paths, 1))
		ratios = np.exp((r-0.5*sigma*sigma)*T/steps + normal_variables * sigma*(T/steps)**0.5)
		# data[:, -1:] select the data block with dims (paths, 1)， while the dims of data[:, -1] is (paths).
		data = np.concatenate((data, data[:, -1:] * ratios), axis=1)

	return data

# def linear_fitting(X, Y):
# 	X = np.array(X)
# 	Y = np.array(Y)
# 	X = np.array([[1]*len(X), X, X*X]).T

# 	coeff = np.linalg.inv(np.matmul(X.T, X))
# 	coeff = np.matmul(coeff, np.matmul(X.T, Y))

# 	Y_fitted = np.matmul(X, coeff)

# 	return Y_fitted

def linear_fitting(X, Y):
	# Linear regression with y = a + bx + cx^2 .
	# Find a, b, c that minimize \sum (y-a-bx-cx^2)^2, there are explict solutions of a, b, c.
	X = np.array(X)
	Y = np.array(Y)
	S0, S1, S2, S3, S4 = len(X), sum(X), sum(X*X), sum(X**3), sum(X**4)
	V0, V1, V2 = sum(Y), sum(Y*X), sum(Y*X*X)

	coeff_mat = np.array([[S0, S1, S2], [S1, S2, S3], [S2, S3, S4]])
	target_vec = np.array([V0, V1, V2])
	inv_coeff_mat = np.linalg.inv(coeff_mat)

	fitting_coeff = np.matmul(inv_coeff_mat, target_vec)

	fitted_Y = fitting_coeff[0]+fitting_coeff[1]*X+fitting_coeff[2]*X*X

	return fitted_Y

def MC_sim_Ame_put_price(r, sigma, S_0, K, T, steps, paths):
	# Sample paths times with "sample_path()" function.
	# Dimensions of data are (paths，steps+1).
	data = sample_paths(r, sigma, S_0, T, steps, paths)

	# Use a list to store the option prices on the ends of each paths.
	option_prices = np.maximum(K-data[:, -1], 0)

	# Calculate the option prices on each path from the end to the front.
	for t_n in range(steps-1, 0, -1):
		# Discount.
		option_prices *= np.exp(-r*T/steps)
		# Re-calculate option price with linear fitting.
		option_prices = linear_fitting(data[:, t_n], option_prices)
		# Decide if we should exercise the option.
		option_prices = np.maximum(option_prices, K-data[:, t_n])
	# Because we can not use linear fitting to the first step in the paths.
	# We need to discount one more time.
	option_prices *= np.exp(-r*T/steps)

	return np.average(option_prices)

if __name__ == "__main__":
	r, sigma, S_0, K, T = 0.1, 0.4, 50, 60, 5/12.
	steps, paths = 150, 20000

	Ame_put_price = MC_sim_Ame_put_price(r, sigma, S_0, K, T, steps, paths)

	print("steps: {}, paths: {} \nAmerican put option price: {:.5f}".format(steps, paths, Ame_put_price))

