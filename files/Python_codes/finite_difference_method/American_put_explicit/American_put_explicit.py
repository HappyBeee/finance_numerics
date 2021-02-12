// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

def American_put_explicit(r, sigma, S_0, K, T, M, N):
	dS = 3*S_0/M
	dt = T/N

	f1 = [max(0.0, K-i*dS) for i in range(M+1)]
	f2 = [None for i in range(M+1)]

	for i in range(N-1, -1, -1):
		f2 = list(f1)
		f1[0] = K
		f1[-1] = 0
		for j in range(1, M, 1):
			f1[j] = f2[j-1]*(-0.5*r*j*dt+0.5*sigma*sigma*j*j*dt)/(1+r*dt)
			f1[j] += f2[j]*(1-sigma*sigma*j*j*dt)/(1+r*dt)
			f1[j] += f2[j+1]*(0.5*r*j*dt+0.5*sigma*sigma*j*j*dt)/(1+r*dt)
		
		for j in range(M+1):
			f1[j] = max(f1[j], K-j*dS)

	
	pos = int(S_0/dS)
	put_price = f1[pos] + (f1[pos+1]-f1[pos])/dS*(S_0-dS*pos)

	return put_price

if __name__ == "__main__":
	put_price = American_put_explicit(0.1, 0.4, 50, 60, 1.0, 200, 20000)
	print("American put price: {0:0.5f}".format(put_price))

