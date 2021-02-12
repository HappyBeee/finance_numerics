// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

def American_put_implicit(r, sigma, S_0, K, T, M, N):
	dS = 3*S_0/M
	dt = T/N

	a = lambda j: 0.5*r*j*dt-0.5*sigma*sigma*j*j*dt
	b = lambda j: 1+r*dt+sigma*sigma*j*j*dt
	c = lambda j: -0.5*r*j*dt-0.5*sigma*sigma*j*j*dt
			 

	f1 = [max(K-i*dS, 0.0) for i in range(M+1)]
	f2 = [None for i in range(M+1)]
	coeffs = np.zeros((M-1, M-1))

	for i in range(N-1, -1, -1):
		f2 = list(f1)
		coeffs[0][0] = b(1)
		coeffs[0][1] = c(1)
		coeffs[M-2][M-2] = b(M-1)
		coeffs[M-2][M-3] = a(M-1) 
		for j in range(2, M-1, 1):
			coeffs[j-1][j-2] = a(j)
			coeffs[j-1][j-1] = b(j)
			coeffs[j-1][j] = c(j)
		coeffs_inv = np.linalg.inv(coeffs)
		Y = f2[1:-1]
		Y[0] -= a(1)*K

		f1[0] = K
		f1[M] = 0
		f1[1:M] = np.matmul(coeffs_inv, Y)

		f1 = np.maximum(f1, K-np.linspace(0, M, M+1)*dS)

	pos = int(S_0/dS)
	put_price = f1[pos] + (f1[pos+1]-f1[pos])/dS*(S_0-pos*dS)

	return put_price

if __name__ == "__main__":
	put_price = American_put_implicit(0.1, 0.4, 50, 50, 5/12.0, 300, 300)
	print("American put price: {0:0.5f}".format(put_price))

