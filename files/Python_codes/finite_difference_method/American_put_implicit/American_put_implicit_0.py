def American_put_implicit(r, sigma, S_0, K, T, M, N):
	dS = 2*S_0/M
	dt = T/N

	f1 = [max(K-i*dS, 0.0) for i in range(M+1)]
	f2 = [None for i in range(M+1)]
	coeffs = [[None, None] for i in range(M+1)]

	for i in range(N-1, -1, -1):
		f2 = list(f1)
		print(coeffs)
		coeffs[0] = [K, 0]
		coeffs[1] = [0, 1]
		for j in range(2, M+1, 1):
			a = 0.5*r*j*dt-0.5*sigma*sigma*j*j*dt
			b = 1+r*dt+sigma*sigma*j*j*dt
			c = -0.5*r*j*dt-0.5*sigma*sigma*j*j*dt
			coeffs[j][0] = f2[j-1]/c - a/c*coeffs[j-2][0] - b/c*coeffs[j-1][0]
			coeffs[j][1] = -a/c*coeffs[j-2][1] - b/c*coeffs[j-1][1]
		x1 = -coeffs[M][0]/coeffs[M][1]
		for j in range(M+1):
			f1[j] = coeffs[j][0] + x1*coeffs[j][1]
			f1[j] = max(K-j*dS, f1[j])
	pos = int(S_0/dS)
	put_price = f1[pos] + (f1[pos+1]-f1[pos])/dS*(S_0-pos*dS)

	return put_price

if __name__ == "__main__":
	put_price = American_put_implicit(0.1, 0.4, 50, 50, 5/12.0, 20, 10)
	print("American put price: {0:0.5f}".format(put_price))

