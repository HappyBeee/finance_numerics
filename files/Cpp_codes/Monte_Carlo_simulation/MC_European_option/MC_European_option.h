// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>

// double normal();

class MC_European_option
{
public:
	MC_European_option(double r, double sigma, double S_0, double K, double T, int N, int M);

	double get_call_price();
	double get_put_price();

	void MC_simulation();
protected:
	double r, sigma, S_0, K, T;
	int N, M;
	double call_price, put_price;
};