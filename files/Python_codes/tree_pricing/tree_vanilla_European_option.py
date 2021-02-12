// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np
import math

E = math.e

class Binomial_tree_sim:
    def __init__(self, r, sigma, S_0, K, T, steps, option_type="european", call_or_put="call"):
        """ Initialize the tree with inputs.
        """
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.K = K
        self.T = T
        self.steps = steps
        
        self.option_type = option_type
        self.call_or_put = call_or_put
        
        # Calculate the tree branching parameters.
        self.dt = self.T/self.steps
        self.u = E**(self.sigma*self.dt**0.5)
        self.d = 1/self.u
        self.p = (E**(self.r*self.dt)-self.d)/(self.u-self.d)
        
        # For storing results.
        self.tree = None
        self.option_price = None
        
        # Build a tree.
        self.build_tree()
        
    def build_tree(self):
        """ Build a tree, which is a list of list.
        """
        self.tree = list()
        for lvl in range(self.steps+1):
            row = list()
            for j in range(lvl+1):
                node = dict()
                node["stock_price"] = self.S_0*self.u**(j)*self.d**(lvl-j)
                node["option_price"] = None
                row.append(node)
            self.tree.append(row)
        return
    
    def calculate_option_price(self):
        """ Calculate option price with given type.
        """
        # Shorter names with local variable.
        r, K, steps = self.r, self.K, self.steps
        dt, p = self.dt, self.p
        
        # Calcualte option prices at the last level on the tree.
        for node in self.tree[-1]:
            # If the option type is "call":
            if self.call_or_put == "call":
                node["option_price"] = max(node["stock_price"]-K, 0)
            # If the option type is "put":
            else:
                node["option_price"] = max(K-node["stock_price"], 0)
        
        # If the option is European type.
        if self.option_type == "european":
            # Calculate backward to get the option price at the root node. 
            for lvl in range(steps-1, -1, -1):
                for j in range(len(self.tree[lvl])):
                    self.tree[lvl][j]["option_price"] = E**(-r*dt)*(p*self.tree[lvl+1][j+1]["option_price"]+\
                                                    (1-p)*self.tree[lvl+1][j]["option_price"])
        
        # If the option is American type, the calculation procedures are similar, 
        # we only need to consider if the option will be exercised.
        else:
            for lvl in range(self.steps-1, -1, -1):
                for j in range(len(self.tree[lvl])):
                    self.tree[lvl][j]["option_price"] = E**(-r*dt)*(p*self.tree[lvl+1][j+1]["option_price"]+\
                                                    (1-p)*self.tree[lvl+1][j]["option_price"])
                    # Decide if the option should be exercised.
                    if self.call_or_put == "call":
                        self.tree[lvl][j]["option_price"] = max(self.tree[lvl][j]["option_price"], \
                                                            self.tree[lvl][j]["stock_price"]-K)
                    else:
                        self.tree[lvl][j]["option_price"] = max(self.tree[lvl][j]["option_price"], \
                                                            K-self.tree[lvl][j]["stock_price"])
        
        self.option_price = self.tree[0][0]["option_price"]

        return 

tree_obj = Binomial_tree_sim(0.05, 0.2, 10, 10, 3, 10, option_type="european", call_or_put="call")
tree_obj.calculate_option_price()
print("European call option price: {:0.4f}".format(tree_obj.option_price))

tree_obj = Binomial_tree_sim(0.05, 0.2, 10, 10, 3, 10, option_type="european", call_or_put="put")
tree_obj.calculate_option_price()
print("European put option price: {:0.4f}".format(tree_obj.option_price))

tree_obj = Binomial_tree_sim(0.05, 0.2, 10, 10, 3, 10, option_type="american", call_or_put="call")
tree_obj.calculate_option_price()
print("American call option price: {:0.4f}".format(tree_obj.option_price))

tree_obj = Binomial_tree_sim(0.05, 0.2, 10, 10, 3, 10, option_type="american", call_or_put="put")
tree_obj.calculate_option_price()
print("American put option price: {:0.4f}".format(tree_obj.option_price))