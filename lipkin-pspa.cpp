#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>

using namespace std;

double epsilon = 1;
int g = 10;
double chi;

double get_epsilon_bar(double sigma0) {return sqrt(pow(epsilon,2)+pow(chi*sigma0,2));}

double get_rpa_freq(double beta, double epsilon_bar, double sigma0)
{
	double omega2 = 4*pow(epsilon_bar,2) - (4*g*epsilon*epsilon*chi/epsilon_bar)*tanh(beta*epsilon_bar/2);
	if(omega2>=0)
		return sqrt(omega2);
	else
	{
		cout << beta << " " << sigma0 << " " << epsilon_bar << endl; return 0.0;
	}
}

double pspa(double x, double beta)
{
	double epsilon_bar = get_epsilon_bar(x);
	double omega = get_rpa_freq(beta,epsilon_bar,x);
	double zeta = pow(2*cosh(beta*epsilon_bar/2), 2*g); 
	double zeta_prime = (omega != 0)? omega/(2*epsilon_bar)*sinh(beta*epsilon_bar)/sinh(beta*omega/2): 2/beta*sinh(beta*epsilon_bar)/(2*epsilon_bar) ;
	
	return exp(-0.5*chi*beta*x*x)*zeta*zeta_prime *sqrt(chi*beta/2*M_PI);
}

double spa(double x, double beta)
{
	double epsilon_bar = get_epsilon_bar(x);
	double zeta = pow(2*cosh(beta*epsilon_bar/2), 2*g); 
	return exp(-0.5*chi*beta*x*x)*zeta * sqrt(chi*beta/2*M_PI);
}

double spa_free_energy(double x, double beta)
{
	double epsilon_bar = get_epsilon_bar(x);
	return 0.5*chi*x*x - 2*g/beta*log(2*cosh(beta*epsilon_bar/2));
}

double sigma_limit = 10;
int gsl_mesh = 10000;

double integrate(double beta, double (*func_x)(double, double))
{
   double up_lim = sigma_limit;
   double low_lim = -sigma_limit;
   double fa, fb,x, step;
   step=(up_lim - low_lim)/((double) gsl_mesh);
   fa=(*func_x)(low_lim, beta);
   fb=(*func_x)(up_lim, beta);
   double trapez_sum=0.0;
   for (int j=1; j <= gsl_mesh-1; j++)
   {
      x=j*step+low_lim;
      trapez_sum+=(*func_x)(x, beta);
   }
   trapez_sum=(trapez_sum+fb+fa)*step;
   return trapez_sum;
}

int main(int argc, char* argv[])
{

	// if(argc!=3) {cerr << "Enter (1) chi (2) beta\n"; exit(1);}
	// chi = atof(argv[1]);
	// double beta = atof(argv[2]);
	// ofstream fout("f_vs_sigma.dat");
	// for(double sigma0 = -10; sigma0 < 10; sigma0 += 0.1)
	// {
	// 	fout << sigma0 << " " << -1/beta*log(spa(sigma0, beta)) << " " << -1/beta*log(pspa(sigma0, beta)) << endl;
	// }

	if(argc!=2) {cerr << "Enter (1) chi \n"; exit(1);}
   chi = atof(argv[1]);
	ofstream fout("lipkin_F_vs_beta.dat");
	for(double beta=0.5; beta <5; beta += 0.1)
	{
		fout << beta << " " << -1/beta*log(integrate(beta,&spa)) << " " << -1/beta*log(integrate(beta, &pspa)) << endl;
	}

	return 0;
}


/* double pspa(double x, void * params)
{
	double beta = *(double *) params;
	return pspa(x,beta);
}
double spa(double x, void * params)
{
	double beta = *(double *) params;
	return spa(x,beta);
} */