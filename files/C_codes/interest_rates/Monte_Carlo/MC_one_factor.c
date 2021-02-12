// All rights reserved, please do not use these codes elsewhere.
// 请不要在任何其它地方使用这些代码，仅供参考。
// 2/10/2021 .

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
        result = rates[p-1][1]+(rates[p][1]-rates[p-1][1])/(rates[p][0]-rates[p-1][0])\
                 *(t-rates[p-1][0]);
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
    for (int i=0; i<steps; i++)
    {
        ti = dt*(1+i);
        theta_t = sigma*sigma/2.0/a*(1.0-exp(-2.0*a*ti));
        theta_t += a*(R0t(rates,nR,ti+dt)*(ti+dt)-R0t(rates,nR,ti)*ti)/dt;
        theta_t += ((ti+dt)*R0t(rates,nR,ti+dt)+(ti-dt)*R0t(rates,nR,ti-dt)\
                   -2.0*ti*R0t(rates,nR,ti))/dt/dt;
        result[i] = theta_t;
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

double * MC_one_factor_bond_option(double K, double L, double T1, double T2, \
                                   double a, double sigma, int steps, int num_paths,\
                                   int intervals, double ** rates, int num_rates)
{
    double dt = T2/(steps+1.0);
    int tp = (int)(T1/dt);
    double dt1 = T1-tp*dt;
    double dt2 = dt-dt1;

    double * theta_ts = Theta_ts(a, sigma, T2, steps, rates, num_rates);
    double R0 = R0t(rates, num_rates, dt);

    double dR = 0.3/intervals;
    double ** Qs = malloc(intervals*sizeof(double *));
    double ** Ps = malloc(intervals*sizeof(double *));
    for (int i=0; i<intervals; i++)
    {
        Qs[i] = malloc(2*sizeof(double));
        Qs[i][0] = 0.0;
        Qs[i][1] = 0.0;
        Ps[i] = malloc(2*sizeof(double));
        Ps[i][0] = 0.0;
        Ps[i][1] = 0.0;
    }

    double * R_path;
    double discount_factor;
    int idx;
    for (int i=0; i<num_paths; i++)
    {
        R_path = MC_one_factor_sample(a, sigma, T2, steps, R0, theta_ts);
        discount_factor = 0.0;
        for (int j=steps; j>tp; j--)
        {
            discount_factor += dt*R_path[j];
        }
        discount_factor += dt2*R_path[tp];
        
        idx = (int) (R_path[tp]/dR);
        idx = (idx<intervals) ? idx : intervals-1;
        idx = (idx>=0) ? idx : 0;
        Ps[idx][0] += 1.0;
        Ps[idx][1] += L*exp(-discount_factor);

        discount_factor = dt1*R_path[tp];
        for (int j=tp-1; j>-1; j--)
        {
            discount_factor += dt*R_path[j];
        }
        Qs[idx][0] += 1.0;
        Qs[idx][1] += exp(-discount_factor);

        free(R_path);
    }

    double bond_price;
    double * call_put_prices = malloc(2*sizeof(double));
    call_put_prices[0] = 0.0;
    call_put_prices[1] = 0.0;
    for(int i=0; i<intervals; i++)
    {
        if (Ps[i][0]<1.E-14)
        {
            continue;
        }
        else
        {
            bond_price = Ps[i][1]/Ps[i][0];
            if (bond_price>K)
            {
                call_put_prices[0] += (bond_price-K)*Qs[i][1]/num_paths;
            }
            else
            {
                call_put_prices[1] += (K-bond_price)*Qs[i][1]/num_paths;
            }
        }
    }

    for (int i=0; i<intervals; i++)
    {
        free(Qs[i]);
        free(Ps[i]);
    }
    free(Qs);
    free(Ps);
    free(theta_ts);

    return call_put_prices;
}



int main()
{
    srand(time(NULL));

    int num_rates = 15;
	double rates_array[15][2] = {{3/365.0, 0.0501722}, {31/365.0, 0.0498284}, {62/365.0, 0.0497234},\
                                {94/365.0, 0.0496157}, {185/365.0, 0.0499058}, {367/365.0, 0.0509389},\
                                {731/365.0, 0.0579733}, {1096/365.0, 0.0630595}, {1461/365.0, 0.0673464},\
                                {1826/365.0, 0.0694816}, {2194/365.0, 0.0708807}, {2558/365.0, 0.0727527}, \
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
    int steps = 500;
    int num_paths = 100000;
    int intervals = 200;

    double * call_put_prices = MC_one_factor_bond_option(K, L, T1, T2, a, sigma, steps,\
                                                         num_paths, intervals, rates, num_rates);

    printf("Steps: %d , number of paths: %d, R(T1) sectors: %d \n", steps, num_paths, intervals);
    printf("Call price: %.6lf , put price: %.6lf \n", call_put_prices[0], call_put_prices[1]);

    for (int i=0; i<num_rates; i++)
    {
    	free(rates[i]);
    }
    free(rates);
    free(call_put_prices);
}