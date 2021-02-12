// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import math
from scipy.stats import norm

# Use the constants in "math" module.
E = math.e
PI = math.pi

# Some parameters as in the example.
r, sigma, S_0, K, T = 0.05, 0.20, 90., 100., 3.
# European option prices function.
def european_option(r, sigma, S_0, K, T):
    """ r: risk free rate.     sigma：asset price volitility.    S_0: asset price.
        K: option exercise price.     T：exercise time of the option.
    """
    
    d_1 = (math.log(S_0/K)+(r+0.5*sigma*sigma)*T)/sigma/T**0.5
    d_2 = d_1 - sigma*T**0.5
    
    call_price = S_0*norm.cdf(d_1) - K*E**(-r*T)*norm.cdf(d_2)
    put_price = K*E**(-r*T)*norm.cdf(-d_2) - S_0*norm.cdf(-d_1)
    
    return (call_price, put_price)

call_price, put_price = european_option(r, sigma, S_0, K, T)
print("European call all price：{0:.5}, European put price{1:.5} 。".format(call_price, put_price))
print("\nput-call parity: \n")
print("   {0:.5} + {1:.5} = {2:.5} + {3:.5}\n".format(call_price, K*E**(-r*T), put_price, S_0))
print("            {0:.5} = {1:.5}\n".format(call_price+K*E**(-r*T), put_price+S_0))