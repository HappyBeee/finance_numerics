import numpy as np

print_interval = 10

def GARCH(data, precision=1.e-3, max_iter=1000, eta=1.e-5):
	N = len(data)
	U2 = [0]*N
	for i in range(1, N):
		U2[i] = (data[i]-data[i-1])/data[i-1]
		U2[i] = U2[i] * U2[i]
	U2[1] += 1.e-8

	variances = [0]*N
	variances[2] = U2[1]
	omega = 1.e-6
	alpha = 0.1
	beta = 0.9

	loss = 0 
	n_iter = 0
	ratio_change = 100
	while n_iter < max_iter and ratio_change > precision:
		n_iter += 1
		loss = 0
		for i in range(2, N, 1):
			variances[i] = omega+alpha*U2[i-1]+beta*variances[i-1]
			# if i%200==0:
			# 	print("~~~~~~~~", i, variances[i])
			loss += np.log(variances[i]) + U2[i]/variances[i]

		d_omega = 0
		d_alpha = 0
		d_beta = 0
		coeff = None
		A = 0
		B = 0
		C = 0
		D = 0
		for i in range(3, N, 1):
			# C, D use A, B from last step.
			C = beta*C + B 
			D = beta*D + A
			# Update A, B.
			A = 1 + beta*A
			B = beta*B + U2[i-1]
			
			# Update omega, alpha, beta changes.
			coeff = 1/variances[i]-U2[i]/variances[i]/variances[i]
			d_omega -= eta/N * coeff * A
			d_alpha -= eta/N * coeff * B
			d_beta -= eta/N * coeff * ((i-2)*beta**(i-3)*variances[2] + alpha*C + omega*D)

			ratio_change = max([abs(d_omega/omega), abs(d_alpha/alpha), abs(d_beta/beta)])

		omega += d_omega*1.e-7
		alpha += d_alpha
		beta += d_beta
		if n_iter%print_interval == 0:
			print("iteration number: ", n_iter, " ,  ", "loss: ", loss, " , parameter ratio change: ", ratio_change)
			print(omega, alpha, beta)

	loss = 0
	for i in range(2, N, 1):
		variances[i] = omega+alpha*U2[i-1]+beta*variances[i-1]
		loss += np.log(variances[i]) + U2[i]/variances[i]
	print("\ntotal iteration times: ", n_iter)
	print("omega, alpha, beta: ", omega, alpha, beta)
	print("total loss: ", loss)

	return (variances, (omega, alpha, beta))




if __name__ == "__main__":
	""" Not stable with this data set, i.e. omega tends to be negative. 
		If set omega to be zero, resulted alpha and beta are similar to EWMA's result.
	"""
	data = np.genfromtxt("GARCHCALCSS&P500.txt", skip_header=1, usecols=(1))
	variances, parameters = GARCH(data, precision=1.e-3, max_iter=20000, eta=1.e-4)
