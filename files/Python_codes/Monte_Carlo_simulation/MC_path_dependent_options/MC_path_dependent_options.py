import numpy as np
import math

E = math.e

def sample_S_path(r, sigma, S_0, T, M):
    S = S_0
    path = [S]
    dt = T/M
    for i in range(M):
        S *= E**((r-0.5*sigma*sigma)*dt+sigma*np.random.normal()*dt**0.5)
        path.append(S)
    return path

def MC_float_lookback_call(r, sigma, S_0, T, M, N):
    call_price = 0
    for i in range(N):
        path = sample_S_path(r, sigma, S_0, T, M)
        call_price += 1.0/N*max(0, path[-1]-min(path))
    return call_price*E**(-r*T)

def MC_average_price_call(r, sigma, S_0, K, T, M, N):
    call_price = 0
    for i in range(N):
        path = sample_S_path(r, sigma, S_0, T, M)
        call_price += 1.0/N*max(0, np.average(path)-K)
    return call_price*E**(-r*T)

if __name__ == "__main__":
    float_lookback_call = MC_float_lookback_call(0.1, 0.4, 50, 1, 100, 10000)
    average_price_call = MC_average_price_call(0.1, 0.4, 50, 50, 1, 100, 10000)
    print("European floating lookback call: {0:.5f}".format(float_lookback_call))
    print("Average price call: {0:.5f}".format(average_price_call))