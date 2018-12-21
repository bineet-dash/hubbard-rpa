#ifndef _ALHASSID_RPA_HPP_INCLUDED
#define _ALHASSID_RPA_HPP_INCLUDED

#include "rpa.hpp"


double rpa_det(double omega, vector <MatrixXd> vt, VectorXd hf, double T)
{
	MatrixXd rpa = MatrixXd::Identity(L,L);
	double mu = get_mu(T, hf);
	for(int alpha=0; alpha<L; alpha++)
	{
		for(int alpha_prime=0; alpha_prime<L; alpha_prime++)
		{
			for(int i=0; i<hf.size(); i++)
			{
				for(int j=0; j<hf.size(); j++)
				{
					rpa(alpha,alpha_prime) += U_prime/2*(vt.at(alpha_prime))(j,i)*(vt.at(alpha))(i,j)*(fermi_fn(hf(i)-mu,T)-fermi_fn(hf(j)-mu,T))/(hf(i)-hf(j)+omega);
				}
			}
		}
	}
	double result = rpa.determinant();
	return result;
}

double modified_bisection(vector<MatrixXd> vt, VectorXd hf, double T, double x_min, double x_max)
{
  int keep_count = 0;
	double phi = (1+sqrt(5))/2;
	double a = x_min;
	double b = x_max;
	double c = b-(b-a)/phi;
	double d = a+ (b-a)/phi;
	
	cout.precision(5);

  while( abs(c-d) > 1e-4)
  {
		// cout << a << " " << b << " " <<  c << " " << d << " \t " << abs(rpa_det(c,vt,hf,T)) << " " << abs(rpa_det(d,vt,hf,T)) <<  endl;
		if(abs(rpa_det(c,vt,hf,T)) < abs(rpa_det(d,vt,hf,T)))
		{
			b = d;
		}
		else
		{
			a = c;
		}
		c = b-(b-a)/phi;
		d = a+ (b-a)/phi;
		keep_count++;
		if(keep_count > 100) break;
	}
	double tentative_root = (a+b)/2.0;
	double f_tent = rpa_det(tentative_root,vt,hf,T);
	if(abs(f_tent) < 0.05)
	{
		return tentative_root;
	}
	else
	{
		cerr << "no root found between " << x_min << " and " << x_max << ". For " << tentative_root << ", Minimum f = " << abs(f_tent) << endl;
		return -1;
	}
}

vector <double> find_poles(vector<MatrixXd> vt, VectorXd hf, double T)
{
	vector <double> poles;
	for(int i=0; i<hf.size(); i++)
	{
		if(i>0 && hf(i)==hf(i-1)) continue;
		for(int j=0; j<i; j++)
		{
			double delta_ij = hf(i)-hf(j);
			double fx = rpa_det(delta_ij+1e-5, vt, hf, T);
			if(fx > 100)
			{
				poles.push_back(delta_ij);
			} 
		}
	}
	sort(poles.begin(), poles.end());

	if(poles.size()>1) //delete duplicates
	{
		auto it=poles.begin();
		while( it!=poles.end()-1)
		{
			if(abs(*(it)-*(it+1)) < 0.01)
			{
				poles.erase(it+1);
			}
			else
			{
				it++;
			}
		}
	}
	return poles;
}

vector <double> find_roots(vector<MatrixXd> vt, VectorXd hf, double T, vector <double>& poles)
{
	vector <double> roots;
	if(poles.size()!=0)
	{
		sort(poles.begin(), poles.end());
		vector <int> irrelevant_poles;
		for(int it=0; it<poles.size(); it++)
		{
			double omega_min = (it==0)?0.0 : poles.at(it-1);
			double omega_max = poles.at(it);
			double root = modified_bisection(vt, hf, T, omega_min, omega_max);
			if(root!= -1)
			{
				roots.push_back(root); 
			} 
			else
			{
				irrelevant_poles.push_back(it);
			}
		}
		for(auto const&i :irrelevant_poles) poles.erase(poles.begin()+i);
	}
	return roots;
}

pair<double, double> rpa_free_energy(VectorXd spa_eivals, vector <double> rpa_roots, vector <double> rpa_poles, double T)
{
  double mu = get_mu(T, spa_eivals);

  double beta = 1/T;
  double spa_part = 0.0, rpa_part_nom = 0.0, rpa_part_denom = 0.0;

  for(int i=0; i<spa_eivals.size(); i++)
  {
    spa_part += (-beta*(spa_eivals(i)-mu) >4.0)?  -beta*(spa_eivals(i)-mu):log(1+exp(-beta*(spa_eivals(i)-mu)));
  }

	if(rpa_roots.size()==rpa_poles.size() && rpa_roots.size() !=0 )
	{
		for(int i=0; i<rpa_roots.size(); i++)
		{
			rpa_part_denom += (rpa_roots.at(i)==0)? log(beta/2): (log_sinh(beta*rpa_roots.at(i)/2) - log(abs(rpa_roots.at(i))) );
		}
		for(int i=0; i<rpa_poles.size(); i++)
		{
			rpa_part_nom += log_sinh(beta*rpa_poles.at(i)/2) - log(abs(rpa_poles.at(i))) ;
		}
	}

  double pspa_F =  -T*(rpa_part_nom- rpa_part_denom)- T*spa_part + mu*L;
	double spa_F = -T*spa_part + mu*L;
	return make_pair(spa_F, pspa_F);
}

pair <double, double> get_free_energy_shortcut(MatrixXd u, VectorXd spa_eivals, double temperature)
{
	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		vt.push_back(u.adjoint()*v_i*u);
	}
	vector <double> poles = find_poles(vt, spa_eivals, temperature);
	vector <double> roots = find_roots(vt, spa_eivals, temperature, poles);
	pair <double, double> free_energies= rpa_free_energy(spa_eivals, roots, poles, temperature);
	return free_energies;
}


VectorXd inttobin(int theValue)
{
  VectorXd v(L);
  for (int i = 0; i < L; ++i)  v(L-1-i) = theValue & (1 << i) ? 1 : 0;
  return v;
}

VectorXd get_field(int i)
{
  VectorXd raw = inttobin(i);
  for(int i=0; i<raw.size(); i++) raw(i) = (raw(i)==0)?-1:1;
  return raw;
}

#endif


/* vector < pair<double,double>> get_fx_list(double temperature, ofstream& fout=nullout)
{
	double omega_low = 0;
	double omega_high = 3*(spa_eivals_global).maxCoeff();
	double basic_omega_step = 0.01;
	double omega_step = basic_omega_step;

	vector <pair<double, double>> fx_list;
	vector <double> Df;// D2fx;

	fx_list.push_back(make_pair(omega_low, rpa_det(omega_low, temperature) ));
	Df.push_back(0);
	fx_list.push_back(make_pair(omega_low, rpa_det(omega_low+omega_step, temperature) ));
	Df.push_back( (fx_list.at(1).second - fx_list.at(0).second)/omega_step );
	double D2fx = abs(Df[1]-Df[0])/omega_step;
	int count_fx=0;

	for(double omega= omega_low+omega_step; omega<omega_high; omega+= omega_step)
	{
		if(abs(D2fx) > 100)
		{
			omega_step = max(1/D2fx, 1e-5);
		}
		else
		{
			omega_step = basic_omega_step;
		}
		count_fx++;

		double f_omega = rpa_det(omega, temperature);
		fx_list.push_back(make_pair(omega,f_omega));
		double Dfx = (fx_list.at(count_fx).second - fx_list.at(count_fx-1).second)/omega_step;
		Df.push_back(Dfx);
		double D2fx = abs(Df[count_fx]-Df[count_fx-1])/omega_step;
		fout << omega << " " << f_omega << " " << Dfx << " " << D2fx << endl;
	}
	return fx_list;
} */

/* double newton_raphson(double (*f)(double, vector <MatrixXd>, VectorXd, double), double T, double min, double max, double min_diff)
{
	auto f_min = f(min,T);
	while (min + min_diff < max)
	{
		auto const mid = 0.5 * min + 0.5 * max;
		auto const f_mid = f(mid,T);

		if ((f_min < 0) == (f_mid < 0)) {
			min = mid;
			f_min = f_mid;
		} else {
			max = mid;
		}
	}
	return min;
} */

// inline double cost(double x){return exp(-abs(x));}

/* double bisection(vector<MatrixXd> vt, VectorXd hf, double T, double x_min, double x_max)
{
  double x_mid, l_mid, r_mid;
  int keep_count = 0;

  while( abs(cost(rpa_det(x_mid,vt,hf,T))-1) > 1e-4)
  {
    x_mid= (x_min+x_max)/2;
    l_mid = (x_min+x_mid)/2;
    r_mid = (x_mid+x_max)/2;
    
		cout << x_min << " " << x_max << " " << x_mid << " " << cost(rpa_det(l_mid,vt,hf,T)) <<  " " << cost(rpa_det(r_mid,vt,hf,T)) << endl;

    if(cost(rpa_det(l_mid,vt,hf,T))> cost(rpa_det(r_mid,vt,hf,T)))
    {
      x_max=x_mid;
    }
    else
    {
      x_min=x_mid;
    }
    keep_count++;
    if(keep_count>100) 
    {
      cerr << "No roots found" << endl;
      return -1;
    }
  }
  return x_mid;
} */
