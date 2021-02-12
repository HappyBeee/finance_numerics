// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import math

E = math.e

class Tree_convertible_bond:
    def __init__(self, r, sigma, S_0, T, lbd, conversion_ratio, callback_price, par_value, recycle_ratio, steps):
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.T = T
        self.lbd = lbd
        self.conversion_ratio = conversion_ratio
        self.callback_price = callback_price
        self.par_value = par_value
        self.recycle_ratio = recycle_ratio
        self.steps = steps
        
        self.dt = self.T/self.steps
        self.u = E**(((self.sigma*self.sigma-self.lbd)*self.dt)**0.5)
        self.d = 1/self.u
        self.p_u = (E**(self.r*self.dt)-self.d*E**(-self.lbd*self.dt))/(self.u-self.d)
        self.p_d = (self.u*E**(-self.lbd*self.dt)-E**(self.r*self.dt))/(self.u-self.d)
        self.p_default = 1-self.p_u-self.p_d
        
        self.bond_price = None
        
        self.tree = None
        
        self.build_tree()
    
    def build_tree(self):
        self.tree = list()
        for lvl in range(self.steps+1):
            row = list()
            for j in range(lvl+1):
                node = dict()
                node["S"] = self.S_0*(self.u**j)*(self.d**(lvl-j))
                node["B"] = None
                row.append(node)
            self.tree.append(row)
        return
    
    def calculate_bond_price(self):
        tree = self.tree
        r, steps = self.r, self.steps
        conversion_ratio, callback_price = self.conversion_ratio, self.callback_price
        recycle_ratio, par_value = self.recycle_ratio, self.par_value
        
        dt, u, d = self.dt, self.u, self.d 
        p_u, p_d, p_default = self.p_u, self.p_d, self.p_default
        
        # Discount factor.
        a = E**(-r*dt)
        
        # Boundary condition.
        for node in tree[-1]:
            node["B"] = max(node["S"]*conversion_ratio, par_value)
        
        # Iteratively calculate back to root node.
        for lvl in range(steps-1, -1, -1):
            for j in range(lvl+1):
                tree[lvl][j]["B"] = a*p_u*tree[lvl+1][j+1]["B"]+a*p_d*tree[lvl+1][j]["B"]
                tree[lvl][j]["B"] += a*p_default*par_value*recycle_ratio
                tree[lvl][j]["B"] = max(min(tree[lvl][j]["B"], callback_price), tree[lvl][j]["S"]*conversion_ratio)
        
        self.bond_price = tree[0][0]["B"]
        
        return
                

tree_obj = Tree_convertible_bond(0.05, 0.3, 50, 0.75, 0.01, 2, 113, 100, 0.4, 10)
tree_obj.calculate_bond_price()
bond_price = tree_obj.bond_price

print("r: 0.05,  sigma: 0.3,  S_0: 50,  T: 0.75,  lambda: 0.01,  conversion_ratio: 2,\n")
print("callback_price: 113,  par_value: 100,  recycle_ratio: 0.4,  steps: 10 . \n\n")
print("convertible bond price: {0:.5f} .".format(bond_price))