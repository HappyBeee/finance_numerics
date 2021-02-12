// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import math

E = math.e

class Tree_lookback_option:
    def __init__(self, r, sigma, S_0, T, steps, strike_price=None, option_type="european"):
        """ Initialize the instance.
        """
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.T = T
        self.steps = steps
        self.option_type = option_type
        
        self.dt = self.T/self.steps
        self.u = E**(self.sigma*self.dt**0.5)
        self.d = 1/self.u
        self.p = (E**(self.r*self.dt)-self.d)/(self.u-self.d)
        
        # If the strike_price is None, the option is floating lookback type, 
        # otherwise it's fixed lookback option (exercise price is strike_price).
        self.strike_price = strike_price
        
        self.call_price = None
        self.put_price = None
        
        # Build a tree.
        self.tree = list()
        self.build_tree()
        
    def build_tree(self):
        """ Build the tree of stock prices. At each node, calculate all possible minimum and maximum stock prices on all
            stock price changing paths that end at this node.
        """
        # Shorten the variable names.
        S_0, steps = self.S_0, self.steps
        u, d, p = self.u, self.d, self.p
        
        self.tree = list()
        for lvl in range(steps+1):
            row = list()
            for j in range(lvl+1):
                node = dict()
                node["S"] = S_0*(u**j)*(d**(lvl-j))
                
                # Stock price historical minimum of maximum with corresponding option prices.
                # Minimum or maximimum price is coded with the power of S_0, d.
                # d is the key in the dict(), value is option price.
                node["Smin_vals"] = dict()
                node["Smax_vals"] = dict()
                
                # Find all the possible hostorical minimums and maximums on each node, 
                # and initialize option prices to be None.
                mins = lvl-j-max(0, lvl-2*j)+1
                maxs = j-max(0, 2*j-lvl)+1
                for k in range(mins):
                    node["Smin_vals"][lvl-j-k] = None
                for k in range(maxs):
                    node["Smax_vals"][-j+k] = None
                
                row.append(node)
            self.tree.append(row)
        
        return 
    
    def calculate_call_price(self):
        """ Calculate lookback call price.
        """
        r, S_0, steps = self.r, self.S_0, self.steps
        dt, u, d, p = self.dt, self.u, self.d, self.p
        strike_price = self.strike_price
        
        # Discount factor for each step.
        a = E**(-r*dt)
        
        # Calculate option prices at the end level of the tree.
        for node in self.tree[-1]:
            if strike_price is None:
                for d_num in node["Smin_vals"]:
                    node["Smin_vals"][d_num] = max(0, node["S"]-S_0*(d**d_num))
            else:
                for d_num in node["Smax_vals"]:
                    node["Smax_vals"][d_num] = max(0, S_0*(d**d_num)-strike_price)
        # Calculate backward to the root node.
        for lvl in range(steps-1, -1, -1):
            for j in range(lvl+1):
                node = self.tree[lvl][j]
                if strike_price is None:
                    for d_num in node["Smin_vals"]:
                        node["Smin_vals"][d_num] = p*a*self.tree[lvl+1][j+1]["Smin_vals"][d_num]
                        node["Smin_vals"][d_num] += (1-p)*a*self.tree[lvl+1][j]["Smin_vals"][max(d_num, lvl-2*j+1)]
                        if self.option_type != "european":
                            node["Smin_vals"][d_num] = max(node["Smin_vals"][d_num], node["S"]-S_0*(d**d_num))
                else:
                    for d_num in node["Smax_vals"]:
                        node["Smax_vals"][d_num] = (1-p)*a*self.tree[lvl+1][j]["Smax_vals"][d_num]
                        node["Smax_vals"][d_num] += p*a*self.tree[lvl+1][j+1]["Smax_vals"][min(d_num, lvl-2*j-1)]
                        if self.option_type != "european":
                            node["Smax_vals"][d_num] = max(node["Smax_vals"][d_num], S_0*(d**d_num)-strike_price)
        
        if strike_price is None:
            self.call_price = self.tree[0][0]["Smin_vals"][0]
        else:
            self.call_price = self.tree[0][0]["Smax_vals"][0]
            
        return 
                    
            
    def calculate_put_price(self):
        """ Calculate lookback put option price.
        """
        r, S_0, steps = self.r, self.S_0, self.steps
        u, d, p = self.u, self.d, self.p
        dt, strike_price = self.dt, self.strike_price
        
        a = E**(-r*dt)
        
        for node in self.tree[-1]:
            if strike_price is None:
                for d_num in node["Smax_vals"]:
                    node["Smax_vals"][d_num] = max(0, S_0*d**d_num-node["S"])
            else:
                for d_num in node["Smin_vals"]:
                    node["Smin_vals"][d_num] = max(0, strike_price-S_0*d**d_num)
        
        for lvl in range(steps-1, -1, -1):
            for j in range(lvl+1):
                node = self.tree[lvl][j]
                if strike_price is None:
                    for d_num in node["Smax_vals"]:
                        node["Smax_vals"][d_num] = (1-p)*a*self.tree[lvl+1][j]["Smax_vals"][d_num]
                        node["Smax_vals"][d_num] += p*a*self.tree[lvl+1][j+1]["Smax_vals"][min(d_num, lvl-2*j-1)]
                        if self.option_type != "european":
                            node["Smax_vals"][d_num] = max(node["Smax_vals"][d_num], S_0*(d**d_num)-node["S"])
                else:
                    for d_num in node["Smin_vals"]:
                        node["Smin_vals"][d_num] = p*a*self.tree[lvl+1][j+1]["Smin_vals"][d_num]
                        node["Smin_vals"][d_num] += (1-p)*a*self.tree[lvl+1][j]["Smin_vals"][max(d_num, lvl-2*j+1)]
                        if self.option_type != "european":
                            node["Smin_vals"][d_num] = max(node["Smin_vals"][d_num], strike_price-S_0*(d**d_num))
        
        if strike_price is None:
            self.put_price = self.tree[0][0]["Smax_vals"][0]
        else:
            self.put_price = self.tree[0][0]["Smin_vals"][0]
        
        return 

print("r = 0.1, sigma = 0.4, S_0 = 50, T = 0.25, steps = 5 .")
tree_obj = Tree_lookback_option(0.1, 0.4, 50, 0.25, 5)
tree_obj.calculate_call_price()
tree_obj.calculate_put_price()
print("European floating lookback call: {0:.5f}".format(tree_obj.call_price))
print("European floating lookback put: {0:.5f}".format(tree_obj.put_price))

tree_obj = Tree_lookback_option(0.1, 0.4, 50, 0.25, 5, option_type="american")
tree_obj.calculate_call_price()
tree_obj.calculate_put_price()
print("American floating lookback call: {0:.5f}".format(tree_obj.call_price))
print("American floating lookback put: {0:.5f}".format(tree_obj.put_price))

print("\nr = 0.1, sigma = 0.4, S_0 = 50, T = 0.25, steps = 5, strike_price = 49 .")
tree_obj = Tree_lookback_option(0.1, 0.4, 50, 0.25, 5, strike_price=49)
tree_obj.calculate_call_price()
tree_obj.calculate_put_price()
print("European fixed lookback call: {0:.5f}".format(tree_obj.call_price))
print("European fixed lookback put: {0:.5f}".format(tree_obj.put_price))

tree_obj = Tree_lookback_option(0.1, 0.4, 50, 0.25, 5, strike_price=49, option_type="american")
tree_obj.calculate_call_price()
tree_obj.calculate_put_price()
print("American fixed lookback call: {0:.5f}".format(tree_obj.call_price))
print("American fixed lookback put: {0:.5f}".format(tree_obj.put_price))