#include "alhassid_rpa.hpp"
#include <chrono>

using namespace std::chrono;

void show_time(milliseconds begin_ms, milliseconds end_ms, string s)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout <<  s << " took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << s << " took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << s << " took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << s << " took " << t << "time. Wrong t received.\n"; }
}

double t=1;
double U_prime=2;
int L=4;

milliseconds begin_ms, end_ms;


int main(int argc, char* argv[])
{
  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) temp, (4) config_int .\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);
	int config_int = atoi(argv[4]);

	vector <MatrixXd> v;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		v.push_back(v_i);
	}

  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());


	MatrixXd sigma = MatrixXd::Zero(L,3);
	sigma.col(2) = get_field(config_int);
  MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  MatrixXd u = spa_spectrum.first.real();
  
  VectorXd spa_eivals = spa_spectrum.second;
	cout << spa_eivals.transpose() << endl << endl;

	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "spa calc ");

	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i_transformed = u.adjoint()*(v.at(it))*u;
		vt.push_back(v_i_transformed);
	}

	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "vt calc ");
	cout << endl;

	cout << sigma.col(2).transpose() << endl;
	cout << rpa_det(0.1,vt, spa_eivals, 0.01) << endl;

	vector <double> poles = find_poles(vt, spa_eivals, temperature);
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "find_poles ");
	cout << endl;

	vector <double> roots = find_roots(vt, spa_eivals, temperature, poles);
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	cout << endl;
	show_time(begin_ms, end_ms, "roots calc ");

	cout << endl << "Poles " << endl;
	for(int it=0; it<poles.size(); it++) cout << poles.at(it) << endl;
  cout << endl;
	cout << "Roots " << endl;
	for(int it=0; it<roots.size(); it++) cout << roots.at(it) << endl;
  cout << endl;
  
	// pair <double, double> free_energies= rpa_free_energy(spa_spectrum.second, roots, poles, temperature);
	// cout << temperature << " " << free_energies.first/L << " " << free_energies.second/L << endl;

	return 0;

}	

/* vector <double> delta;
for(int i=0; i<spa_eivals.size(); i++)
{
	for(int j=0; j<i; j++)
	{
		double delta_ij = spa_eivals(i)-spa_eivals(j);
		delta.push_back(delta_ij);
	}
	sort(delta.begin(), delta.end());
}
{
	auto it=delta .begin();
	while( it!=delta.end()-1)
	{
		if(abs(*(it)-*(it+1)) < 0.01)	delta.erase(it+1);
		else it++;
	}
	for(const auto& i:delta) cout << i << "\n"; cout << endl;
} */