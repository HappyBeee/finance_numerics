// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

def get_loss(U2, omega, alpha, beta):
	sigma2 = U2[1]
	loss = 0
	for i in range(2, len(U2)):
		loss += np.log(sigma2)+U2[i]/sigma2
		sigma2 = omega+beta*sigma2+alpha*U2[i]
	return loss

def GARCH_optimal_parameters(data, M1=80, M2=30):
	# Two steps exhaustive search.
	# First step set omega = VL(1-alpha-beta).
	# Second step, probe around the optimal omega, alpha, beta from step one.
	N = len(data)-1
	U2 = [0]*(1+N)
	for i in range(1, N+1):
		U2[i] = (data[i]-data[i-1])/data[i-1]
		U2[i] = U2[i]*U2[i]
	U2 = np.array(U2)
	VL = np.average(U2[1:])

	min_loss = float("inf")
	loss = None
	opt1_omega = None
	opt1_alpha = None
	opt1_beta = None
	for i in range(M1):
		beta = 0.5 + i*0.5/M1
		for j in range(M1):
			alpha = 0.01 + j*0.5/M1
			if alpha + beta >= 1:
				continue
			omega = VL*(1-alpha-beta)
			loss = get_loss(U2, omega, alpha, beta)
			if loss < min_loss:
				min_loss = loss
				opt1_omega = omega
				opt1_alpha = alpha
				opt1_beta = beta
	print("VL = ", VL)
	print("Optimal alpha, beta = ", opt1_alpha, opt1_beta)
	print("Omega = VL(1-alpha-beta) = ", opt1_omega)
	print("Total loss: ", min_loss)
	print("\n\n")

	min_loss = float("inf")
	loss = None
	opt2_omega = None
	opt2_alpha = None
	opt2_beta = None
	for i in range(M2):
		beta = opt1_beta-0.025+i*0.05/M2
		for j in range(M2):
			alpha = opt1_alpha-0.025+j*0.05/M2
			if alpha+beta >= 1:
				continue
			for k in range(M2):
				omega = 0.5*opt1_omega+k*opt1_omega/M2
				loss = get_loss(U2, omega, alpha, beta)
				if loss < min_loss:
					min_loss = loss
					opt2_omega = omega
					opt2_alpha = alpha
					opt2_beta = beta
	print("Final omega, alpha, beta = ", opt2_omega, opt2_alpha, opt2_beta)
	print("Total loss: ", min_loss)
	print("\n\n")

	return (opt2_omega, opt2_alpha, opt2_beta)

def GARCH_predict(data, omega, alpha, beta):
	N = len(data)-1 # Data index ranges from 0 to N.
	U2 = [0]*(1+N)
	for i in range(1, N+1):
		U2[i] = (data[i]-data[i-1])/data[i-1]
		U2[i] = U2[i]*U2[i]

	variances = [0, 0, U2[1]]
	for i in range(2, N+1):
		variances.append(omega+beta*variances[-1]+alpha*U2[i])

	return variances


if __name__ == "__main__":
	data = np.genfromtxt("GARCHCALCSS&P500.txt", skip_header=1, usecols=(1))

	omega, alpha, beta = GARCH_optimal_parameters(data, 70, 20)
