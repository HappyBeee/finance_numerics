// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import math

E = math.e

class Tree_asian_option:
    def __init__(self, r, sigma, S_0, K, T, steps, points):
        """ Initialize the instance, "points" is the number of historical average stock prices
            at each nodes.
        """
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.K = K
        self.T = T
        self.steps = steps
        self.points = points
        
        self.dt = self.T/self.steps
        self.u = E**(self.sigma*self.dt**0.5)
        self.d = 1/self.u
        self.p = (E**(self.r*self.dt)-self.d)/(self.u-self.d)
        
        self.call_price = None
        self.tree = list()
        self.build_tree()
    
    def get_max_min_ave(self, i, j):
        """ Calculate the minimum and maximum historical average stock prices at node (i,j).
        """
        u, d, S_0 = self.u, self.d, self.S_0
        max_ave = S_0*(1-u**(j+1))/(1-u)+S_0*(u**j)*d*(1-d**(i-j))/(1-d)
        max_ave /= 1+i
        min_ave = S_0*(1-d**(i-j+1))/(1-d)+S_0*(d**(i-j))*u*(1-u**j)/(1-u)
        min_ave /= 1+i
        
        return (max_ave, min_ave)
    
    def interpolation(self, x, ref_list):
        """ Linear interpolation, "ref_list" a number list with dimension (points, 2), 
            "ref_list" is defaulted to be sorted as an increasing list by values in {ref_list[i][0]} .
        """
        left, right = 0, len(ref_list)-1
        pos = 0
        
        # If out of range, return the boundary values.
        # In fact, this would not happen on the tree we built. 
        if x > ref_list[right][0]:
            return ref_list[right][1]
        if x < ref_list[left][0]:
            return ref_list[left][1]
        
        # Find the position of x in {ref_list[i][0]} by binary search.
        while left < right:
            pos = int((left+right)/2)
            if x == ref_list[pos][0]:
                return ref_list[pos][1]
            if x > ref_list[pos][0]:
                left = pos+1
            else:
                right = pos-1
        if x > ref_list[left][0]:
            pos = left+1
        else:
            pos = left
        
        # Linear interpolation.
        result = (ref_list[pos][1]-ref_list[pos-1][1])/(ref_list[pos][0]-ref_list[pos-1][0])
        result *= (x-ref_list[pos-1][0])
        result += ref_list[pos-1][1]
        
        return result
    
    def build_tree(self):
        S_0, steps, points = self.S_0, self.steps, self.points
        u, d, p = self.u, self.d, self.p
        
        self.tree = list()
        for lvl in range(steps+1):
            row = list()
            for j in range(lvl+1):
                node = dict()
                node["S"] = S_0*(u**j)*(d**(lvl-j))
                node["F_S"] = list()
                max_ave, min_ave = self.get_max_min_ave(lvl, j)
                for k in range(points):
                    node["F_S"].append([(max_ave-min_ave)/(points-1)*k+min_ave, None])
                row.append(node)
            self.tree.append(row)
        
        return
    
    def calculate_call_price(self):
        """ Calculate the European average price call.
        """
        r, S_0, K = self.r, self.S_0, self.K
        steps, points = self.steps, self.points
        dt, u, d, p = self.dt, self.u, self.d, self.p
        
        a = E**(-r*dt)
        # Boundary condition.
        for node in self.tree[-1]:
            for k in range(points):
                node["F_S"][k][1] = max(0, node["F_S"][k][0]-K)
        
        # Calculate backward to the root node.
        for lvl in range(steps-1, -1, -1):
            for j in range(lvl+1):
                node = self.tree[lvl][j]
                for k in range(points):
                    new_ave_u = (node["F_S"][k][0]*(lvl+1)+self.tree[lvl+1][j+1]["S"])/(lvl+2)
                    new_ave_d = (node["F_S"][k][0]*(lvl+1)+self.tree[lvl+1][j]["S"])/(lvl+2)
                    option_price_u = self.interpolation(new_ave_u, self.tree[lvl+1][j+1]["F_S"])
                    option_price_d = self.interpolation(new_ave_d, self.tree[lvl+1][j]["F_S"])
                    node["F_S"][k][1] = a*p*option_price_u+a*(1-p)*option_price_d
        
        self.call_price = self.tree[0][0]["F_S"][0][1]
        
        return    

print("r = 0, sigma = 0.4, S_0 = 50, K = 50, T = 1, steps = 60, points = 100 .\n")
tree_obj = Tree_asian_option(0.1, 0.4, 50, 50, 1, 60, 100)
tree_obj.calculate_call_price()
print("European average price call: {0:.5f}".format(tree_obj.call_price))