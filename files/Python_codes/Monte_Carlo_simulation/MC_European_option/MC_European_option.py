import numpy as np
import math 

class Monte_Carlo_European_option:
    def __init__(self, r, sigma, S_0, K, T, M, N):
        self.r = r
        self.sigma = sigma
        self.S_0 = S_0
        self.K = K
        self.T = T
        self.M = M
        self.N = N
        self.call_price = None
        self.put_price = None
    
    def MC_simulation(self):
        self.call_price = 0
        self.put_price = 0
        dt = self.T/self.M
        
        # Sample N stock price changing paths.
        for i in range(self.N):
            S = self.S_0
            # M samples of Normal distribution for the stock price path.
            for j in range(self.M):
                S = S*(1+self.r*dt+self.sigma*dt**0.5*np.random.normal())
                # Another option for discretization.
                # S = S*math.e**((self.r-0.5*self.sigma*self.sigma)*dt+self.sigma*dt**0.5*np.random.normal())
            self.call_price += max(0, S-self.K)/self.N
            self.put_price += max(0, self.K-S)/self.N
        
        self.call_price *= math.e**(-self.r*self.T)
        self.put_price *= math.e**(-self.r*self.T)
        
        return

if __name__ == "__main__":
    MC_obj = Monte_Carlo_European_option(0.05, 0.2, 90, 100, 1, 400, 40000)
    MC_obj.MC_simulation()
    print("European call option price：{0:.5f}".format(MC_obj.call_price))
    print("European put option price：{0:.5f}".format(MC_obj.put_price))