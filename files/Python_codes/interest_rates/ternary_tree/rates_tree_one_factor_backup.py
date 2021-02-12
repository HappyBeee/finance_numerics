import math

init_rates = [[3/365, 0.0501722], [31/365, 0.0498284], [62/365, 0.0497234], [94/365, 0.0496157],\
              [185/365, 0.0499058], [367/365, 0.0509389], [731/365, 0.0579733], [1096/365, 0.0630595], \
              [1461/365, 0.0673464], [1826/365, 0.0694816], [2194/365, 0.0708807], [2558/365, 0.0727527], \
              [2922/365, 0.0730852], [3287/365, 0.0739790], [3653/365, 0.0749015]]

def calculate_init_zero_bonds(init_rates, dt, steps):
	""" Calculate the initial zero bonds prices for calibrating the rates tree.
		Use linear interpolation for zero rate at each time.
		Time range form ti = 0, to ti = (steps+1)*dt. N-steps branching has N+1 time intervals.
	"""
	init_rates = list(init_rates)
	P_ti = [1]
	if len(init_rates) == 1:
		r_flat = init_rates[0][1]
		for i in range(1, steps+2):
			P_ti.append(math.exp(-r_flat*i*dt))

		return P_ti

	# Extend init rates to t=0.
	r_0 = init_rates[0][1]-(init_rates[1][1]-init_rates[0][1])/(init_rates[1][0]-init_rates[0][0])\
			*init_rates[0][0]
	r_0 = max(r_0, 0)
	init_rates = [[0, r_0]] + init_rates
	pos = 0
	for i in range(1, steps+2):
		ti = i*dt
		while ti > init_rates[pos][0]:
			pos += 1
		r_i = init_rates[pos-1][1] + (init_rates[pos][1]-init_rates[pos-1][1])\
			  /(init_rates[pos][0]-init_rates[pos-1][0])*(ti-init_rates[pos-1][0])
		P_ti.append(math.exp(-r_i*ti))

	return P_ti

def rates_tree(a, sigma, T, steps, init_rates):
	""" Ternary tree based on Hull White one factor model.
		N steps branching has N+1 time intervals.
		Initial rates is of form [[t_i, r_i], ...], has time range greater than T.
	"""
	dt = T/(steps+1)
	dR = sigma*(3*dt)**0.5
	# Max radius for the Ternary tree according to John Hull's book.
	max_radius = int(math.ceil(0.184/a/dt))
	# Zero bonds prices from initial rates.
	P_ti = calculate_init_zero_bonds(init_rates, dt, steps)

	# Initialize the tree as a list of lists of nodes.
	# Node is a dictionary.
	tree = [[{"Q": 1, "Probs":[None, None, None], "prob": 1, "rate": None}]]
	# Build tree from root level, by
	#   1) Construct next level's empty nodes. 
	#	2) Calculate current level's interest rates, branching probabilities.
	#   3) Sending current level nodes' Qs to next level.
	for lvl in range(steps):
		next_lvl_nodes = [{"Q":0, "Probs":[None, None, None], "prob": 0, "rate": None}\
						 for _ in range(2*min(max_radius, lvl+1)+1)]
		self.tree.append(next_lvl_nodes)

		nodes = tree[lvl]
		# Calibrate interest rate by matching initial bonds prices.
		alpha = math.log(sum([nodes[j]["Q"]*math.exp(-j*dR*dt) for j in range(-min(lvl, max_radius), min(lvl, max_radius)+1)]))
		alpha = (alpha - math.log(P_ti[lvl+1]))/dt
		for j in range(-min(lvl, max_radius-1), min(lvl, max_radius-1)+1):
			nodes[j]["rate"] = j*dR+alpha
			nodes[j]["Probs"][0] = 1./6+0.5*(a*a*j*j*dt*dt+a*j*dt)
			nodes[j]["Probs"][1] = 2./3-a*a*j*j*dt*dt
			nodes[j]["Probs"][2] = 1./6+0.5*(a*a*j*j*dt*dt-a*j*dt)
			for k in range(3):
				tree[lvl+1][j+k-1]["Q"] += nodes[j]["Probs"][k]*nodes[j]["Q"]*math.exp(-nodes[j]["rate"]*dt)
				tree[lvl+1][j+k-1]["prob"] += nodes[j]["Probs"][k]*nodes[j]["prob"]
		if lvl >= max_radius:
			for sign in [-1, 1]:
				j_b = max_radius*sign
				nodes[j_b]["rate"] = j_b*dR+alpha
				nodes[j_b]["Probs"][1+sign] = 7./6+0.5*(a*a*j_b*j_b*dt*dt-3*sign*a*j_b*dt)
				nodes[j_b]["Probs"][1] = -1./3-a*a*j_b*j_b*dt*dt+2*sign*a*j_b*dt
				nodes[j_b]["Probs"][1-sign] = 1./6+0.5*(a*a*j_b*j_b*dt*dt-sign*a*j_b*dt)
				for k in range(3):
					tree[lvl+1][j_b+k-sign-1]["Q"] += nodes[j_b]["Probs"][k]*nodes[j_b]["Q"]*math.exp(-nodes[j_b]["rate"]*dt)
					tree[lvl+1][j_b+k-sign-1]["prob"] += nodes[j_b]["Probs"][k]*nodes[j_b]["prob"]
	
	lvl += 1
	alpha = math.log(sum([tree[lvl][j]["Q"]*math.exp(-j*dR*dt) for j in range(-min(lvl, max_radius), min(lvl, max_radius)+1)]))
	alpha = (alpha - math.log(P_ti[lvl+1]))/dt
	for j in range(-min(lvl, max_radius), mmin(lvl, max_radius)+1):
		tree[lvl][j]["rate"] = j*dR+alpha

	return tree