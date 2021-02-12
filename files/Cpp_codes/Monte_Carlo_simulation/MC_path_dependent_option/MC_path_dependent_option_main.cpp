// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

#include "MC_path_dependent_option.cpp"

using namespace std;

int main(){
	srand(time(nullptr));

	double r = 0.1, sigma = 0.4;
	double S_0 = 50, K = 50, T = 1;
	int N = 200, M = 40000;

	double float_lookback_call_price = float_lookback_call(r, sigma, S_0, T, N, M);
	double average_price_call_price = average_price_call(r, sigma, S_0, K, T, N, M);
	cout << "Float lookback call option price: " << float_lookback_call_price << endl;
	cout << "Average price call option price: " << average_price_call_price << endl;
}