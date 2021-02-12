#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double normal()
{
    double x, y;
    double result;
    x = (0.5+rand())/((double)RAND_MAX+1.0);
    y = (0.5+rand())/((double)RAND_MAX+1.0);
    result = sqrt(-2.0*log(x))*cos(6.283185307179586*y);
    return result;
}

double R0t(double ** rates, int num_rates, double t)
{
    double result;
    if (t<=rates[0][0])
    {
        result = rates[0][1];
    }
    else if (t>=rates[num_rates-1][0])
    {
        result = rates[num_rates-1][1];
    }
    else
    {
        int p = 0;
        while (rates[p][0]<t)
        {
            p++;
        }
        result = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])*(t-rates[p-1][0]);
    }
    return result;
}

double * Theta_ts(double a, double sigma, double T, int steps, double ** rates, int num_rates)
{
    int nR = num_rates;
    double dt = T/(1.0+steps);
    double * result = malloc(steps*sizeof(double));
    double theta_t;
    double ti;
    for (int i=1; i<steps+1; i++)
    {
        ti = dt*i;
        theta_t = sigma*sigma/2.0/a*(1.0-exp(-2.0*a*ti));
        theta_t += a*(R0t(rates,nR,ti+dt)*(ti+dt)-R0t(rates,nR,ti)*ti)/dt;
        theta_t += ((ti+dt)*R0t(rates,nR,ti+dt)+(ti-dt)*R0t(rates,nR,ti-dt)-2.0*ti*R0t(rates,nR,ti))/dt/dt;
        result[i-1] = theta_t;
    }
    return result;
}

double * MC_one_factor_sample(double a, double sigma, double T, int steps, double R0, double * theta_ts)
{
    double dt = T/(1.0+steps);
    double * result = malloc((steps+1)*sizeof(double));
    result[0] = R0;
    double Ri = R0;
    for (int i=1; i<steps+1; i++)
    {
        Ri = Ri+(theta_ts[i-1]-a*Ri)*dt+sigma*sqrt(dt)*normal();
        result[i] = Ri;
    }
    return result;
}

double * MC_one_factor_bond_option(double K, double L, double t, double T, double a, double sigma,\
                                   int steps, int num_paths, int intervals, double ** rates, int num_rates)
{
    double dt = T/(steps+1.0);
    int tp = (int) (t/dt);
    double dt1 = t-tp*dt;
    double dt2 = dt-dt1;

    double * theta_ts = Theta_ts(a, sigma, T, steps, rates, num_rates);
    double R0 = R0t(rates, num_rates, dt);

    double dR = 0.3/intervals;
    double ** discounts = malloc(intervals*sizeof(double *));
    double ** bond_prices = malloc(intervals*sizeof(double *));
    for (int i=0; i<intervals; i++)
    {
        discounts[i] = malloc(2*sizeof(double));
        discounts[i][0] = 0.0;
        discounts[i][1] = 0.0;
        bond_prices[i] = malloc(2*sizeof(double));
        bond_prices[i][0] = 0.0;
        bond_prices[i][1] = 0.0;
    }

    double * R_path;
    double discount_factor;
    int interval_idx;
    for (int i=0; i<num_paths; i++)
    {
        R_path = MC_one_factor_sample(a, sigma, T, steps, R0, theta_ts);
        discount_factor = 0.0;
        for (int j=steps; j>tp; j--)
        {
            discount_factor += dt*R_path[j];
        }
        discount_factor += dt2*R_path[tp];
        
        interval_idx = (int) (R_path[tp]/dR);
        interval_idx = (interval_idx<intervals) ? interval_idx : intervals-1;
        interval_idx = (interval_idx>=0) ? interval_idx : 0;
        bond_prices[interval_idx][0] += 1.0;
        bond_prices[interval_idx][1] += L * exp(-discount_factor);

        discount_factor = dt1*R_path[tp];
        for (int j=tp-1; j>-1; j--)
        {
            discount_factor += dt*R_path[j];
        }
        discounts[interval_idx][0] += 1.0;
        discounts[interval_idx][1] += exp(-discount_factor);

        free(R_path);
    }

    double bond_price;
    double * call_put_prices = malloc(2*sizeof(double));
    call_put_prices[0] = 0.0;
    call_put_prices[1] = 0.0;
    for(int i=0; i<intervals; i++)
    {
        if (bond_prices[i][0]<1.E-14)
        {
            continue;
        }
        else
        {
            bond_price = bond_prices[i][1]/bond_prices[i][0];
            if (bond_price>K)
            {
                call_put_prices[0] += (bond_price-K)*discounts[i][1]/num_paths;
            }
            else
            {
                call_put_prices[1] += (K-bond_price)*discounts[i][1]/num_paths;
            }
        }
    }

    // for (int i=0; i<intervals; i++)
    // {
    //     bond_price = (bond_prices[i][0]>1.E-14) ? bond_prices[i][1]/bond_prices[i][0] : 0.0;
    //     double discount = (discounts[i][0]>1.E-14) ? discounts[i][1]/discounts[i][0] : 0.0;
    //     printf("Rt: %.6lf  counts: %.0lf  bond_price: %.5lf  discount: %.5lf \n", \
    //             i*dR, bond_prices[i][0], bond_price, discount);
    // }
    // printf("\n\n");


    for (int i=0; i<intervals; i++)
    {
        free(discounts[i]);
        free(bond_prices[i]);
    }
    free(discounts);
    free(bond_prices);
    free(theta_ts);

    return call_put_prices;
}



int main()
{
    srand(time(NULL));

	int num_rates = 15;
	double rates_array[15][2] = {{3/365.0, 0.0501722}, {31/365.0, 0.0498284}, {62/365.0, 0.0497234}, {94/365.0, 0.0496157},\
              			  {185/365.0, 0.0499058}, {367/365.0, 0.0509389}, {731/365.0, 0.0579733}, {1096/365.0, 0.0630595}, \
              			  {1461/365.0, 0.0673464}, {1826/365.0, 0.0694816}, {2194/365.0, 0.0708807}, {2558/365.0, 0.0727527}, \
              			  {2922/365.0, 0.0730852}, {3287/365.0, 0.0739790}, {3653/365.0, 0.0749015}};
    double ** rates = malloc(num_rates*sizeof(double *));
    for (int i=0; i<num_rates; i++)
    {
    	rates[i] = malloc(2*sizeof(double));
    	rates[i][0] = rates_array[i][0];
    	rates[i][1] = rates_array[i][1];
    }

    double a = 0.1;
    double sigma = 0.01;
    double T1 = 3.0;
    double T2 = 9.0;
    double K = 63.0;
    double L = 100.0;
    int steps = 1000;
    int num_paths = 100000;
    int intervals = 2000;

    double * call_put_prices = MC_one_factor_bond_option(K, L, T1, T2, a, sigma, steps, num_paths, intervals, rates, num_rates);

    printf("Steps: %d , number of paths: %d, Rt sectors: %d \n", steps, num_paths, intervals);
    printf("Call price: %.6lf , put price: %.6lf \n", call_put_prices[0], call_put_prices[1]);

    for (int i=0; i<num_rates; i++)
    {
    	free(rates[i]);
    }
    free(rates);
    free(call_put_prices);
}