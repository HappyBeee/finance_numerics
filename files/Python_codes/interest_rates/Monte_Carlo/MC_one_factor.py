// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

INIT_RATES = [[3/365, 0.0501722], [31/365, 0.0498284], [62/365, 0.0497234], [94/365, 0.0496157],\
              [185/365, 0.0499058], [367/365, 0.0509389], [731/365, 0.0579733], [1096/365, 0.0630595], \
              [1461/365, 0.0673464], [1826/365, 0.0694816], [2194/365, 0.0708807], [2558/365, 0.0727527], \
              [2922/365, 0.0730852], [3287/365, 0.0739790], [3653/365, 0.0749015]]


def R0t(init_rates, t):
	rates = list(init_rates)
	result = None
	if t <= rates[0][0]:
		result = rates[0][1]
	elif t >= rates[-1][0]:
		result = rates[-1][1]
	else:
		p = 0
		while rates[p][0] < t:
			p += 1
		result = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])*(t-rates[p-1][0])
	return result

def Theta_ts(a, sigma, T, steps, init_rates):
	rates = list(init_rates)
	dt = T/(1.0+steps)
	result = []
	theta_t = None
	for i in range(1, steps+1, 1):
		ti = dt*i
		theta_t = sigma*sigma/2.0/a*(1.0-np.exp(-2.0*a*ti))
		theta_t += a*(R0t(rates, ti+dt)*(ti+dt)-R0t(rates, ti)*ti)/dt
		theta_t += ((ti+dt)*R0t(rates, ti+dt)+(ti-dt)*R0t(rates, ti-dt)-2.0*ti*R0t(rates, ti))/dt/dt
		result.append(theta_t)
	return result

def MC_one_factor_sample(a, sigma, T, steps, R0, theta_ts):
	dt = T/(1.0+steps)
	result = [R0]
	Ri = R0
	for i in range(1, steps+1, 1):
		Ri = Ri+(theta_ts[i-1]-a*Ri)*dt+sigma*np.sqrt(dt)*np.random.normal()
		result.append(Ri)
	return result

def MC_bond_option_one_factor(K, L, T1, T2, a, sigma, steps, num_paths, intervals, init_rates):
	rates = list(init_rates)
	dt = T2/(1.0+steps)
	tp = int(T1/dt)
	dt1 = T1-tp*dt
	dt2 = dt-dt1

	theta_ts = Theta_ts(a, sigma, T2, steps, rates)
	R0 = R0t(rates, dt)

	dR = 0.3/intervals
	discounts = [[] for _ in range(intervals)]
	bond_prices = [[] for _ in range(intervals)]

	R_path = None
	discount_factor = None
	Rts_idx = None
	for i in range(num_paths):
		R_path = MC_one_factor_sample(a, sigma, T2, steps, R0, theta_ts)
		discount_factor = 0.0
		for j in range(steps, tp, -1):
			discount_factor += dt*R_path[j]
		discount_factor += dt2*R_path[tp]

		Rts_idx = int(R_path[tp]/dR)
		Rts_idx = min(intervals-1, Rts_idx)
		Rts_idx = max(0, Rts_idx)
		bond_prices[Rts_idx].append(L*np.exp(-discount_factor))

		discount_factor = dt1*R_path[tp]
		for j in range(tp-1, -1, -1):
			discount_factor += dt*R_path[j]
		discounts[Rts_idx].append(np.exp(-discount_factor))

	bond_price = None
	call_price, put_price = 0.0, 0.0
	for i in range(intervals):
		if len(bond_prices[i]) == 0:
			continue
		else:
			bond_price = np.average(bond_prices[i])
			if bond_price > K:
				call_price += (bond_price-K)*np.sum(discounts[i])/num_paths
			else:
				put_price += (K-bond_price)*np.sum(discounts[i])/num_paths

	return call_price, put_price 


if __name__ == "__main__":
	a, sigma = 0.1, 0.01
	K, L = 63.0, 100.0
	T1, T2 = 3.0, 9.0
	steps = 200
	num_paths = 20000
	intervals = 200

	call_price, put_price = MC_bond_option_one_factor(K, L, T1, T2, a, sigma, steps, num_paths, intervals, INIT_RATES)
	print("\n\ncall, put prices: {0:.5f}  {1:.5f}".format(call_price, put_price))


