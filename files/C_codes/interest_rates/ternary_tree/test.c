#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Node
{
	int lvl;
	int j;
	double prob;
	double Probs[3];
	// Probabilities of branching, Probs[0]:P_d, Probs[1]:P_m, Probs[2]:P_u.
	double R;
	double Q;
};


struct Rates_tree
{
	int steps;
	double T;
	double a;
	double sigma;
	double ** init_rates;
	
	struct Node ** tree;
};

int main()
{
	int num_init_rates = 15;
	double rates[15][2] = {{3/365.0, 0.0501722}, {31/365.0, 0.0498284}, {62/365.0, 0.0497234}, {94/365.0, 0.0496157},\
              			  {185/365.0, 0.0499058}, {367/365.0, 0.0509389}, {731/365.0, 0.0579733}, {1096/365.0, 0.0630595}, \
              			  {1461/365.0, 0.0673464}, {1826/365.0, 0.0694816}, {2194/365.0, 0.0708807}, {2558/365.0, 0.0727527}, \
              			  {2922/365.0, 0.0730852}, {3287/365.0, 0.0739790}, {365.03/365.0, 0.0749015}};
    double ** init_rates = malloc(15*sizeof(double*));
    for (int i=0; i<num_init_rates; i++)
    {
    	init_rates[i] = malloc(2*sizeof(double));
    	init_rates[i][0] = rates[i][0];
    	init_rates[i][1] = rates[i][1];
    }

}