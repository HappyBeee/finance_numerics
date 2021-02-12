// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

import numpy as np

def EWMA(values, precision=1.e-3):
    M = int(1/precision)
    values = np.array(values)

    U = (values[1:]-values[:-1])/values[:-1]
    U_squared = U*U
    
    opt_lbd = None
    min_loss = float("inf")

    for i in range(1, M):
        lbd = float(i)/M
        sigma_squared = U_squared[0]
        loss = 0
        for j in range(1, len(U_squared)):
            loss += np.log(sigma_squared)+U_squared[j]/sigma_squared
            sigma_squared = lbd * sigma_squared + (1-lbd)*U_squared[j]
        if loss < min_loss:
            min_loss = loss
            opt_lbd = lbd
    
    print("Optimal lambda: ", opt_lbd)
    print("Total loss: ", min_loss)
    Vars = [0, U_squared[0]]
    for i in range(1, len(U_squared)):
        Vars.append(Vars[-1]*opt_lbd+(1-opt_lbd)*U_squared[i])
    
    return (Vars, opt_lbd)

if __name__ == "__main__":
    data = np.genfromtxt("GARCHCALCSS&P500.txt", skip_header=1, usecols=(1))
    #data = np.genfromtxt("EURUSDExchangerates.txt", skip_header=1, usecols=(1))
    Vars, lbd = EWMA(data, precision=1.e-3)