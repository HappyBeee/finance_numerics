// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

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

double * calculate_init_zero_bonds(double ** init_rates, double T, int steps)
{
	double dt = T/(steps+1.0);
	double * Pti = malloc((steps+2) * sizeof(double));
	int p = 0;
	double ti, ri;
	for (int i=0; i<steps+2; i++)
	{
		ti = dt*i;
		while (ti>init_rates[p][0])
		{
			p++;
		}
		if (p==0)
		{
			ri = init_rates[0][1];
		}
		else
		{
			ri = init_rates[p-1][1]+(init_rates[p][1]-init_rates[p-1][1])\
				 /(init_rates[p][0]-init_rates[p-1][0])*(ti-init_rates[p-1][0]);
		}
		Pti[i] = exp(-ri*ti);;
	}
	return Pti;
}


struct Rates_tree * build_rates_tree(double a, double sigma, double T, int steps, double ** init_rates, int num_init_rates)
{
	struct Rates_tree * rates_tree = malloc(sizeof(struct Rates_tree)); 
	// Store rates tree parameters.
	rates_tree->a = a;
	rates_tree->sigma = sigma;
	rates_tree->T = T;
	rates_tree->steps = steps;
	rates_tree->tree = malloc((steps+1)*sizeof(struct Node *));

	rates_tree->init_rates = malloc(num_init_rates*sizeof(double *));
	for (int i=0; i<num_init_rates; i++)
	{	
		rates_tree->init_rates[i] = malloc(2*sizeof(double));
		rates_tree->init_rates[i][0] = init_rates[i][0];
		rates_tree->init_rates[i][1] = init_rates[i][1];
	}

	// Parameters for building tree.
	double dt = T/(steps+1.0);
	double dR = sigma*sqrt(3.0*dt);
	int j_max = (int)ceil(0.184/a/dt);
	double * Pti = calculate_init_zero_bonds(rates_tree->init_rates, T, steps);

	// Shorten variable name.
	struct Node ** tree = rates_tree->tree;
	
	tree[0] = malloc(1*sizeof(struct Node));
	tree[0][0].lvl = 0;
	tree[0][0].j = 0;
	tree[0][0].prob = 1.0;
	tree[0][0].Q = 1.0;
	tree[0][0].R = -1.E100;
	tree[0][0].Probs[0] = -1.E100;
	tree[0][0].Probs[1] = -1.E100;
	tree[0][0].Probs[2] = -1.E100;

	int N1 = steps<j_max ? steps : j_max; // Number of the first standard branching levels. 

	// A Node type pointer to replace long node names.
	struct Node * node;

	double alpha, S;
	// First build tree before nonstandard branching.
	for (int lvl=0; lvl<N1; lvl++)
	{
		// Build next level's empty nodes.
		tree[lvl+1] = malloc((2*lvl+3)*sizeof(struct Node));
		for (int j=-lvl-1; j<lvl+2; j++)
		{
			node = &tree[lvl+1][j+lvl+1];
			node->lvl = lvl+1;
			node->j = j;
			node->prob = 0;
			node->Q = 0;
			node->R = -1.E100;
			node->Probs[0] = -1.E100;
			node->Probs[1] = -1.E100;
			node->Probs[2] = -1.E100;
		}

		// Process with current level's nodes.
		//  1. Calculate R = dR*j+alpha.
		//  2. Calculate branching probabilities.
		//  3. Send probabilities of arriving at nodes and Qs to next level.
		S = 0;
		for (int j=-lvl; j<lvl+1; j++)
		{
			S += tree[lvl][j+lvl].Q * exp(-j*dR*dt);
		}
		alpha = (log(S)-log(Pti[lvl+1]))/dt;

		for (int j=-lvl; j<lvl+1; j++)
		{
			node = &tree[lvl][j+lvl];
			node->R = j*dR+alpha;
			node->Probs[0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt);
			node->Probs[1] = 2.0/3.0-a*a*j*j*dt*dt;
			node->Probs[2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt);
			
			for (int k=0; k<3; k++)
			{
				tree[lvl+1][j+lvl+k].prob += node->prob * node->Probs[k];
				tree[lvl+1][j+lvl+k].Q += node->Q * node->Probs[k] * exp(-node->R*dt);
			}
		}	
	}

	// Consider if the tree has nonstandard branching.
	// Process with standard nodes first, then nonstandard nodes.
	if (steps<=j_max)
	{}
	else
	{
		for (int lvl=j_max; lvl<steps; lvl++)
		{
			tree[lvl+1] = malloc((2*j_max+1)*sizeof(struct Node));
			for (int j=-j_max; j<j_max+1; j++)
			{
				node = &tree[lvl+1][j+j_max];
				node->lvl = lvl+1;
				node->j = j;
				node->prob = 0.0;
				node->Q = 0.0;
				node->R = -1.E100;
				node->Probs[0] = -1.E100;
				node->Probs[1] = -1.E100;
				node->Probs[2] = -1.E100;
			}

			S = 0;
			for (int j=-j_max; j<j_max+1; j++)
			{
				S += tree[lvl][j+j_max].Q*exp(-j*dR*dt);
			}
			alpha = (log(S)-log(Pti[lvl+1]))/dt;

			// Deal with standard nodes first.
			for (int j=-j_max+1; j<j_max; j++)
			{
				node = &tree[lvl][j+j_max];

				node->R = j*dR+alpha;
				node->Probs[0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt);
				node->Probs[1] = 2.0/3.0-a*a*j*j*dt*dt;
				node->Probs[2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt);

				for (int k=0; k<3; k++)
				{
					tree[lvl+1][j+j_max+k-1].prob += node->prob * node->Probs[k];
					tree[lvl+1][j+j_max+k-1].Q += node->Q * node->Probs[k] * exp(-node->R*dt);
				}
			}
			// Deal with bottom and top nonstandard nodes.
			int j = -j_max;
			node = &tree[lvl][0];
			node->R = j*dR+alpha;
			node->Probs[0] = 7.0/6.0+0.5*(a*a*j*j*dt*dt+3.0*a*j*dt);
			node->Probs[1] = -1.0/3.0-a*a*j*j*dt*dt-2.0*a*j*dt;
			node->Probs[2] = 1.0/6.0+0.5*(a*a*j*j*dt*dt+a*j*dt);
			for (int k=0; k<3; k++)
			{
				tree[lvl+1][k].prob += node->prob * node->Probs[k];
				tree[lvl+1][k].Q += node->Q * node->Probs[k] * exp(-node->R*dt);
			}
			j = j_max;
			node = &tree[lvl][2*j_max];
			node->R = j*dR+alpha;
			node->Probs[0] = 1.0/6.0+0.5*(a*a*j*j*dt*dt-a*j*dt);
			node->Probs[1] = -1.0/3.0-a*a*j*j*dt*dt+2.0*a*j*dt;
			node->Probs[2] = 7.0/6.0+0.5*(a*a*j*j*dt*dt-3.0*a*j*dt);
			for (int k=0; k<3; k++)
			{
				tree[lvl+1][2*j_max-2+k].prob += node->prob * node->Probs[k];
				tree[lvl+1][2*j_max-2+k].Q += node->Q * node->Probs[k] * exp(-node->R*dt);
			}
		}
	}

	// Fill last level's rates.
	S = 0;
	int radius = steps<j_max ? steps : j_max;
	for (int j=-radius; j<radius+1; j++)
	{
		S += tree[steps][j+radius].Q * exp(-j*dR*dt);
	}
	alpha = (log(S)-log(Pti[steps+1]))/dt;
	for (int j=-radius; j<radius+1; j++)
	{
		tree[steps][j+radius].R = j*dR+alpha;
	}

	free(Pti);

	return rates_tree;
}

void free_rates_tree(struct Rates_tree * rates_tree, int num_init_rates)
{
	for (int i=0; i<rates_tree->steps+1; i++)
	{
		free(rates_tree->tree[i]);
	}
	for (int i=0; i<num_init_rates; i++)
	{
		free(rates_tree->init_rates[i]);
	}
	free(rates_tree->tree);
	free(rates_tree->init_rates);
	free(rates_tree);
}

double P(double t, double ** init_rates)
{
	double r;
	int p = 0;
	while (t>init_rates[p][0])
	{
		p++;
	}
	if (p==0)
	{
		r = init_rates[0][1];
	}
	else
	{
		r = init_rates[p-1][1]+(init_rates[p][1]-init_rates[p-1][1])\
			/(init_rates[p][0]-init_rates[p-1][0])*(t-init_rates[p-1][0]);
	}
	return exp(-r*t);
}

double B(double t, double T, double a)
{
	return (1.0-exp(-a*(T-t)))/a;
}

double PtT_explicit(double t, double T, double R, double dt, double a, double sigma, double ** init_rates)
{
	double A_hat = P(T, init_rates)/P(t, init_rates);
	A_hat *= exp(-log(P(t+dt, init_rates)/P(t, init_rates))*B(t, T, a)/B(t, t+dt, a));
	A_hat /= exp(sigma*sigma/4.0/a*(1.0-exp(-2.0*a*t))*B(t, T, a)*(B(t, T, a)-B(t, t+dt, a)));

	return A_hat * exp(-B(t, T, a)/B(t, t+dt, a)*dt*R);
}

double * bond_option(double K, double L, double t, double T, int steps, double a, 
					 double sigma, double ** init_rates, int num_init_rates)
{
	double * call_put_prices = malloc(2*sizeof(double));
	call_put_prices[0] = 0.0;
	call_put_prices[1] = 0.0;

	double dt = t/steps;
	struct Rates_tree * rates_tree = build_rates_tree(a, sigma, t+dt, steps, init_rates, num_init_rates);
	struct Node ** tree = rates_tree->tree;

	double bond;
	int width = 1-2*tree[steps][0].j;
	for (int j=0; j<width; j++)
	{
		bond = L * PtT_explicit(t, T, tree[steps][j].R, dt, a, sigma, init_rates);
		if (bond>K)
		{
			call_put_prices[0] += (bond-K) * tree[steps][j].Q;
		}
		else
		{
			call_put_prices[1] += (K-bond) * tree[steps][j].Q;
		}
	}

	free_rates_tree(rates_tree, num_init_rates);

	return call_put_prices;
}

void print_tree(struct Rates_tree * rates_tree)
{
	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	for (int i=0; i<rates_tree->steps+1; i++)
	{
		printf("\ni = %d\n R: ", i);
		for (int j=0; j<1-2*rates_tree->tree[i][0].j; j++)
		{
			printf("%.4lf  ", rates_tree->tree[i][j].R);
		}
		printf("\n Q: ");
		for (int j=0; j<1-2*rates_tree->tree[i][0].j; j++)
		{
			printf("%.4lf  ", rates_tree->tree[i][j].Q);
		}
	}
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

int main()
{
	int num_init_rates = 15;
	double rates[15][2] = {{3/365.0, 0.0501722}, {31/365.0, 0.0498284}, {62/365.0, 0.0497234},\
							{94/365.0, 0.0496157}, {185/365.0, 0.0499058}, {367/365.0, 0.0509389},\
							{731/365.0, 0.0579733}, {1096/365.0, 0.0630595}, {1461/365.0, 0.0673464},\
							{1826/365.0, 0.0694816}, {2194/365.0, 0.0708807}, {2558/365.0, 0.0727527}, \
              			  {2922/365.0, 0.0730852}, {3287/365.0, 0.0739790}, {3653/365.0, 0.0749015}};
    double ** init_rates = malloc(15*sizeof(double*));
    for (int i=0; i<num_init_rates; i++)
    {
    	init_rates[i] = malloc(2*sizeof(double));
    	init_rates[i][0] = rates[i][0];
    	init_rates[i][1] = rates[i][1];
    }

    double a = 0.1;
    double sigma = 0.01;
    double t = 3.0;
    double T = 9.0;
    double K = 63.0;
    double L = 100.0;
    int steps = 200;

    double * call_put_prices = bond_option(K, L, t, T, steps, a, sigma, init_rates, num_init_rates);
    
    printf("Interest rate tree branching steps: 200 .\n");
    printf("Call price: %.5f , Put price: %.5f .\n", call_put_prices[0], call_put_prices[1]);

    // steps = 15;
    // struct Rates_tree * rates_tree = build_rates_tree(a, sigma, t, steps, init_rates, num_init_rates);
    // print_tree(rates_tree);
    // free_rates_tree(rates_tree, num_init_rates);

    for (int i=0; i<num_init_rates; i++)
    {
    	free(init_rates[i]);
    }
    free(init_rates);
    free(call_put_prices);

    // printf("size of int: %ld \n", sizeof(int));
    // printf("size of double: %ld \n", sizeof(double));
    // printf("size of int * : %ld \n", sizeof(int *));
    // printf("size of double ** : %ld \n", sizeof(double **));
    // printf("size of Node: %ld \n", sizeof(struct Node));
    // printf("size of Rates_tree: %ld \n", sizeof(struct Rates_tree));
    // printf("size of Node *: %ld \n", sizeof(struct Node *));
    // printf("size of Rates_tree **: %ld \n", sizeof(struct Rates_tree **));

    // Size of struct increases as a multiplication of 8, 8 * (ceil(num(int)/2.0)+nums(double)) .
    // Size of pointer type variable has size 8.
}