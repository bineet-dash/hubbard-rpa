#include <functional>
#include <cmath>
#include <iostream>
#include "rpa.hpp"


double t=1;
double U_prime=2;
int L=4;

double find_root(std::function<double(double)> f, double min, double max, double epsilon)
{
	auto f_min = f(min);
	while (min + epsilon < max)
	{
		auto const mid = 0.5 * min + 0.5 * max;
		auto const f_mid = f(mid);

		if ((f_min < 0) == (f_mid < 0)) {
			min = mid;
			f_min = f_mid;
		} else {
			max = mid;
		}
	}
	return min;
}

double rpa_det(vector <MatrixXd> v, VectorXd hf, double omega, double T)
{
	MatrixXd rpa = MatrixXd::Identity(L,L);
	
	for(int alpha=0; alpha<L; alpha++)
	{
		for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
		{
			for(int i=0; i<hf.size(); i++)
			{
				for(int j=0; j<hf.size(); j++)
				{			
					rpa(alpha,alpha_prime) += U_prime/2*(v.at(alpha_prime))(j,i)*(v.at(alpha))(i,j)*(fermi_fn(hf(i),T)-fermi_fn(hf(j),T))/(hf(i)-hf(j)+omega);
				}
			}

		}
	}

	double result = rpa.determinant();
	return result;
}

vector <double> find_poles(vector <pair<double,double>> fx_list)
{
  vector <double> poles;
  bool WATCH= false;
  double x_left, x_right;
  double F_MAX_BOUND = 1e+3;
  
	for(auto it=fx_list.begin(); it != fx_list.end(); it++)
	{
		double x = (*it).first; double fx= (*it).second;

		if(abs(fx) > F_MAX_BOUND && !WATCH)
		{
			x_left = x;
			WATCH =true;
		}

		if(WATCH && abs(fx) < F_MAX_BOUND) 
		{
			x_right = x;
			WATCH = false;
			double pole = (x_left+x_right)/2.0;
			poles.push_back(pole); cout << pole << endl;
		}
	}
  return poles;
}

vector <double> find_roots(vector <pair<double,double>> fx_list)
{
  vector <double> roots;
  bool WATCH= false;
  double x_left, x_right;
  double F_MIN_BOUND = 1e-2;

	for(auto it=fx_list.begin(); it != fx_list.end(); it++)
	{
		double x = (*it).first; double fx= (*it).second;
    
    if(abs(fx) < F_MIN_BOUND)
		{
			WATCH = true;
		}
    else 
    {
      WATCH = false;
    }

    if(WATCH) 
    {
      double fx_next = (*(it+1)).second;
      double fx_prev = (*(it-1)).second;
      if( fx*fx_next < 0)
      {
        roots.push_back(x);
        // WATCH = false;
        cout << x << endl;
      } 
      else if((fx-fx_prev)*(fx_next-fx) < 0)
      {
        roots.push_back(x);
        // WATCH = false;
        cout << x << endl;
      }
    }
  }

}

double rpa_free_energy(VectorXd spa_eivals, vector <double> rpa_eivals, vector <double> poles, double T)
{
  double mu = get_mu(T, spa_eivals);
  spa_eivals = spa_eivals.array()-mu;

  double beta = 1/T;
  double spa_part = 0.0, rpa_part_nom = 0.0, rpa_part_denom = 0.0;

  for(int i=0; i<spa_eivals.size(); i++) 
  {
    spa_part += (-beta*spa_eivals(i)>4.0)?  -beta*spa_eivals(i):log(1+exp(-beta*spa_eivals(i)));
  }

  for(int i=0; i<rpa_eivals.size(); i++)
  {
    rpa_part_denom += (rpa_eivals.at(i)==0)? log(beta/2): (log_sinh(beta*rpa_eivals.at(i)/2) - log(abs(rpa_eivals.at(i))) );
  }
  
  for(int i=0; i<poles.size(); i++)
  {
    rpa_part_nom += (poles.at(i)==0)? log(beta/2): (log_sinh(beta*poles.at(i)/2) - log(abs(poles.at(i))) );
  }

  return -T*(2*rpa_part_nom-2*rpa_part_denom+spa_part)/L + mu;
}

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U, (3) temp.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);

	MatrixXd sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  // for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
  MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  MatrixXd u = spa_spectrum.first.real();

	vector <MatrixXd> v;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		v.push_back(v_i);
	}

	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i_transformed = u.adjoint()*(v.at(it))*u;
		vt.push_back(v_i_transformed);
	}
  
	double omega_low = 0;//2*(spa_spectrum.second).minCoeff();
	double omega_high = 2*(spa_spectrum.second).maxCoeff();
	double omega_step = 0.01;

	ofstream fout("omega.dat");

	vector <pair<double, double>> fx_list;
	for(double omega= omega_low; omega<omega_high; omega+= omega_step)
	{	
		double f_omega = rpa_det(vt, spa_spectrum.second, omega, temperature);
		fx_list.push_back(make_pair(omega,f_omega));	
		fout << omega << " " << f_omega << endl;
	}

  vector <double> poles = find_poles(fx_list);
  vector <double> roots = find_roots(fx_list);



}	