// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#define E 2.718281828459045
#define PI 3.141592653589793

struct Node
{
	double asset_price;
	double up_prob;
	double down_prob;

	double node_call_price;
	double node_put_price;
};


class Tree_vanilla_option
{
public:
	Tree_vanilla_option(double r, double sigma, double S_0, double K, double T, int steps);
	
	// user end functions
	Node ** get_tree();
	double get_call_price();
	double get_put_price();

protected:
	// necessary parameters
	int steps;
	double r, sigma, S_0, K, T;
	
	// main targets
	Node ** tree;

	void build_tree();
	void calculate_call_price();
	void calculate_put_price();
};


void print_tree(Node **, int, bool);

