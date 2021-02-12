import numpy as np

INIT_RATES = [[3/365, 0.0501722], [31/365, 0.0498284], [62/365, 0.0497234], [94/365, 0.0496157],\
              [185/365, 0.0499058], [367/365, 0.0509389], [731/365, 0.0579733], [1096/365, 0.0630595], \
              [1461/365, 0.0673464], [1826/365, 0.0694816], [2194/365, 0.0708807], [2558/365, 0.0727527], \
              [2922/365, 0.0730852], [3287/365, 0.0739790], [3653/365, 0.0749015]]

def R_0t(init_rates, t):
	rates = list(init_rates)
	p = 0
	R = None
	while rates[p][0] < t:
		p += 1
	if p == 0:
		R = rates[0][1]
	else:
		R = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])*(t-rates[p-1][0])
	return R

def get_theta_ts(a, sigma, T, steps, init_rates):
	rates = list(init_rates)
	dt = T/(1.0+steps)
	theta_ts = []
	theta_t = None
	for i in range(1, steps+1, 1):
		ti = dt*i
		theta_t = sigma*sigma/2.0/a*(1.0-np.exp(-2.0*ti*a))
		theta_t += a/dt*(R_0t(rates, ti+dt)*(ti+dt)-R_0t(rates, ti)*ti)
		theta_t += 1.0/dt/dt*((ti+dt)*R_0t(rates, ti+dt)+(ti-dt)*R_0t(rates, ti-dt)-2.0*ti*R_0t(rates, ti))
		theta_ts.append(theta_t)
	return theta_ts

def MC_one_factor_sample(a, sigma, T, steps, theta_ts, R0):
	dt = T/(1.0+steps)
	Rs = [R0]
	Ri = R0
	for i in range(1, steps+1, 1):
		Ri = Ri+(theta_ts[i-1]-a*Ri)*dt+sigma*np.sqrt(dt)*np.random.normal()
		Rs.append(Ri)
	return Rs

def MC_one_factor_bond_option(K, L, t, T, a, sigma, steps, num_paths, init_rates):
	call_price = 0
	put_price = 0
	dt = T/(1.0+steps)
	rates = list(init_rates)
	tp = int(t/dt)
	dt1 = t-tp*dt
	dt2 = dt-dt1

	R0 = R_0t(rates, dt)
	theta_ts = get_theta_ts(a, sigma, T, steps, rates)
	bond_price = None
	option_price = None
	discount_factor = None
	for i in range(num_paths):
		bond_price = L
		option_price = None
		discount_factor = 0
		R_path = MC_one_factor_sample(a, sigma, T, steps, theta_ts, R0)
		for j in range(steps, tp, -1):
			discount_factor += dt*R_path[j]
		discount_factor += dt2*R_path[tp]
		bond_price /= np.exp(discount_factor)

		option_price = bond_price - K

		discount_factor = dt1*R_path[tp]
		for j in range(tp-1, -1, -1):
			discount_factor += dt*R_path[j]
		option_price /= np.exp(discount_factor)
		if option_price > 0:
			call_price += option_price
		else:
			put_price -= option_price

	return call_price/float(num_paths), put_price/float(num_paths)

if __name__ == "__main__":
	K, L = 63.0, 100.0
	t, T = 3.0, 9.0
	a, sigma = 0.1, 0.01
	steps = 1000
	num_paths = 2000

	call_price, put_price = MC_one_factor_bond_option(K, L, t, T, a, sigma, steps, num_paths, INIT_RATES)
	print("Bond option call and put prices: {0:.5f}  {1:.5f}".format(call_price, put_price))