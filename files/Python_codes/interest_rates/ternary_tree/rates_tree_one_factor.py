// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

init_rates = [[3/365, 0.0501722], [31/365, 0.0498284], [62/365, 0.0497234], [94/365, 0.0496157],\
              [185/365, 0.0499058], [367/365, 0.0509389], [731/365, 0.0579733], [1096/365, 0.0630595], \
              [1461/365, 0.0673464], [1826/365, 0.0694816], [2194/365, 0.0708807], [2558/365, 0.0727527], \
              [2922/365, 0.0730852], [3287/365, 0.0739790], [3653/365, 0.0749015]]


def calculate_init_zero_bonds(init_zero_rates, T, steps):
	""" Expect T is smaller than init_zero_rates time range.
	"""
	rates = list(init_zero_rates)
	dt = T/(steps+1)
	Pti = [1]
	if len(rates) == 1:
		r_flat = rates[0][1]
		for i in range(1, steps+2, 1):
			Pti.append(np.exp(-r_flat*i*dt))

	# Add a r_0 at time=0, by extending rates curve to left.
	r_0 = rates[0][1]-(rates[1][1]-rates[0][1])/(rates[1][0]-rates[0][0])*rates[0][0]
	rates = [[0, r_0]] + rates

	p = 0
	for i in range(1, steps+2, 1):
		ti = dt*i
		while ti > rates[p][0]:
			p += 1
		ri = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])*(ti-rates[p-1][0])
		Pti.append(np.exp(-ri*ti))

	return Pti


def build_rates_tree(a, sigma, T, steps, init_zero_rates):
	dt = T/(steps+1.0)
	dR = sigma*np.sqrt(3*dt)
	j_max = int(np.ceil(0.184/a/dt))

	Pti = calculate_init_zero_bonds(init_zero_rates, T, steps)
	tree = [[{"Probs":[None, None, None], "prob":1, "Q":1, "R":None}]]
	
	# First build tree levels before nonstandard branching.
	for lvl in range(min(steps, j_max)):
		# Build next level's empty nodes.
		next_lvl_nodes = []
		for j in range(2*lvl+3):
			next_lvl_nodes.append({"Probs":[None, None, None], "prob":0, "Q":0, "R":None})
		tree.append(next_lvl_nodes)

		# Process with current level's nodes.
		#  1. Calculate R = dR*j+alpha.
		#  2. Calculate branching probabilities.
		#  3. Send probabilities of arriving at nodes and Qs to next level.
		nodes = tree[lvl]
		S = 0
		for j in range(-lvl, lvl+1, 1):
			S += nodes[j]["Q"]*np.exp(-j*dR*dt)
		alpha = (np.log(S)-np.log(Pti[lvl+1]))/dt

		for j in range(-lvl, lvl+1, 1):
			nodes[j]["R"] = j*dR+alpha

			nodes[j]["Probs"][0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt)
			nodes[j]["Probs"][1] = 2.0/3.0-a*a*j*j*dt*dt
			nodes[j]["Probs"][2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt)

			for k in range(3):
				tree[lvl+1][j+k-1]["prob"] += nodes[j]["prob"]*nodes[j]["Probs"][k]
				tree[lvl+1][j+k-1]["Q"] += nodes[j]["Q"]*nodes[j]["Probs"][k]*np.exp(-nodes[j]["R"]*dt)

	# Consider if the tree has nonstandard branching.
	# Process with standard nodes first, then nonstandard nodes.
	if steps > j_max:
		for lvl in range(j_max, steps, 1):
			next_lvl_nodes = []
			for j in range(2*j_max+1):
				next_lvl_nodes.append({"Probs":[None, None, None], "prob":0, "Q":0, "R":None})
			tree.append(next_lvl_nodes)

			nodes = tree[lvl]
			S = 0
			for j in range(-j_max, j_max+1, 1):
				S += nodes[j]["Q"]*np.exp(-j*dR*dt)
			alpha = (np.log(S)-np.log(Pti[lvl+1]))/dt

			# Deal with standard nodes first.
			for j in range(-j_max+1, j_max, 1):
				nodes[j]["R"] = j*dR+alpha

				nodes[j]["Probs"][0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt)
				nodes[j]["Probs"][1] = 2.0/3.0-a*a*j*j*dt*dt
				nodes[j]["Probs"][2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt)

				for k in range(3):
					tree[lvl+1][j+k-1]["prob"] += nodes[j]["prob"]*nodes[j]["Probs"][k]
					tree[lvl+1][j+k-1]["Q"] += nodes[j]["Q"]*nodes[j]["Probs"][k]*np.exp(-nodes[j]["R"]*dt)

			# Deal with bottom and top nonstandard nodes.
			j = -j_max
			nodes[j]["R"] = j*dR+alpha
			nodes[j]["Probs"][0] = 7.0/6.0+0.5*(a*a*j*j*dt*dt+3*a*j*dt)
			nodes[j]["Probs"][1] = -1.0/3.0-a*a*j*j*dt*dt-2.0*a*j*dt
			nodes[j]["Probs"][2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt)
			for k in range(3):
				tree[lvl+1][j+k]["prob"] += nodes[j]["prob"]*nodes[j]["Probs"][k]
				tree[lvl+1][j+k]["Q"] += nodes[j]["Q"]*nodes[j]["Probs"][k]*np.exp(-nodes[j]["R"]*dt)
			j = j_max
			nodes[j]["R"] = j*dR+alpha
			nodes[j]["Probs"][0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt) 
			nodes[j]["Probs"][1] = -1.0/3.0-a*a*j*j*dt*dt+2.0*a*j*dt
			nodes[j]["Probs"][2] = 7.0/6.0+0.5*(a*a*j*j*dt*dt-3*a*j*dt)
			for k in range(3):
				tree[lvl+1][j+k-2]["prob"] += nodes[j]["prob"]*nodes[j]["Probs"][k]
				tree[lvl+1][j+k-2]["Q"] += nodes[j]["Q"]*nodes[j]["Probs"][k]*np.exp(-nodes[j]["R"]*dt)
	
	# Fill last level's rates.
	S = 0
	radius = min(steps, j_max)
	for j in range(-radius, radius+1, 1):
		S += tree[steps][j]["Q"]*np.exp(-j*dR*dt)
	alpha = (np.log(S)-np.log(Pti[steps+1]))/dt
	for j in range(-radius, radius+1, 1):
		tree[steps][j]["R"] = j*dR+alpha

	return tree


def PtT_explicit(t, T, R, dt, a, sigma, init_zero_rates):
	rates = list(init_zero_rates)
	def r(t, rates):
		# Find the interest rate (yield) for zero bond.
		p = 0
		while t > rates[p][0]:
			p += 1
		if p == 0:
			r_t = rates[0][1]
		else:
			r_t = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])*(t-rates[p-1][0])
		return r_t

	def P(t, rates):
		return np.exp(-r(t, rates)*t)

	def B(t, T):
		return (1.0-np.exp(-a*(T-t)))/a

	B_hat = B(t, T)/B(t, t+dt)*dt
	A_hat = P(T, rates)/P(t, rates)/np.power(P(t+dt, rates)/P(t, rates), B(t, T)/B(t, t+dt))
	A_hat /= np.exp(sigma*sigma/4/a*(1-np.exp(-2*a*t))*B(t, T)*(B(t, T)-B(t, t+dt)))

	return A_hat * np.exp(-B_hat*R)


def future_bond_option(K, L, t, T, steps, a, sigma, init_zero_rates):
	""" Build a rates tree from t=0 to t=t+dt.
		Use P(t, T) explicit expression for Hull-White one factor model at 
		 the end of the tree.
	""" 
	dt = t/steps
	tree = build_rates_tree(a, sigma, t+dt, steps, init_zero_rates)

	call_price = 0
	put_price = 0
	for j in range(len(tree[-1])):
		P = L*PtT_explicit(t, T, tree[-1][j]["R"], dt, a, sigma, init_zero_rates)
		if P > K:
			call_price += (P-K)*tree[-1][j]["Q"]
		else:
			put_price += (K-P)*tree[-1][j]["Q"]

	return call_price, put_price


if __name__ == "__main__":
	call_price, put_price = future_bond_option(63, 100, 3, 9, 500, 0.1, 0.01, init_rates)
	print("Future bond option, call and put prices: ", call_price, put_price)

	tree = build_rates_tree(0.1, 0.01, 9, 500, init_rates)

	# tree = build_rates_tree(0.1, 0.01, 3, 15, init_rates)
	# for i in range(16):
	# 	print("i = ", i)
	# 	rs = []
	# 	radius = (len(tree[i])-1)//2
	# 	for j in range(-radius, radius+1):
	# 		rs.append([tree[i][j]["R"], tree[i][j]["Q"]])
	# 	print(rs)
