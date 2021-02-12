// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

#include <iostream>
#include <cmath>

#define E 2.71828182846

struct Node
{
	double S;

	// The array of pairs of stock average price and related option price. 
	double ** F_S;
};

class Tree_asian_option
{
public:
	Tree_asian_option(double r, double sigma, double S_0, double K, double T, int steps, int points);

	Node ** get_tree();
	double get_call_price();
	double interpolation(double x, double ** ref_list);

protected:
	double r, sigma, S_0, K, T;
	int steps, points;

	double dt, u, d, p;

	Node ** tree;
	void build_tree();
	void calculate_call_price();

	// Get linear interpolation value.
	// Get max and min averages of historical path stock prices. 
	double * get_mean_range(int lvl, int j);
};
