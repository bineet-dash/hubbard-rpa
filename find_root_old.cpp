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

ofstream nullout;

vector < pair<double,double>> get_fx_list(double temperature, vector <MatrixXd> vt, VectorXd hf, ofstream& fout=nullout)
{
	double omega_low = 0;//2*(spa_spectrum.second).minCoeff();
	double omega_high = 3.0*hf.maxCoeff();
	double OMEGA_STEP = 0.01;
	vector <pair<double, double>> fx_list;
	for(double omega= omega_low; omega<omega_high; omega+= OMEGA_STEP)
	{
		double f_omega = rpa_det(omega, vt, hf, temperature);
		fx_list.push_back(make_pair(omega,f_omega));
		fout << omega << " " << f_omega << endl;
		cout << omega << " done. \r"; cout.flush();
	}
	return fx_list;
}


vector <double> find_roots_old(vector <pair<double,double>> fx_list, double T)
{
  double F_MIN_BOUND = 0.5;
  double TOLERANCE = 1e-2;
  vector <double> roots;

	for(auto it=fx_list.begin(); it != fx_list.end(); it++)
	{
    double fx= (*it).second;
    double fx_next = (*(it+1)).second;
    double fx_prev = (*(it-1)).second;

    if(abs(fx) < F_MIN_BOUND)
		{
      if( fx*fx_next < 0)
      {
        // double bisected_root= newton_raphson(&rpa_det, T, (*it).first, (*(it+1)).first, 1e-11);
        double bisected_root = ((*it).first+(*(it+1)).first)/2.0;
        roots.push_back(bisected_root);
      }

			if((fx-fx_prev)*(fx_next-fx) < 0 && abs(fx)<TOLERANCE)
			{
				roots.push_back((*it).first);
			}
    }
  }

	return roots;
}

vector <double> find_poles_old(vector <pair<double,double>> fx_list)
{
  bool WATCH_POLE = false;
  double F_MAX_BOUND = 10*fx_list.at(1).second;
	double OMEGA_STEP = 0.01;
  
	vector <double> poles;
  double x_left, x_right;

	for(auto it=fx_list.begin(); it != fx_list.end(); it++)
	{
		double x = (*it).first; double fx= (*it).second;

		if(abs(fx) > F_MAX_BOUND && !WATCH_POLE)
		{
			x_left = x;
			WATCH_POLE =true;
		}

		if(WATCH_POLE && abs(fx) < F_MAX_BOUND)
		{
			x_right = x;
			WATCH_POLE = false;
			double pole = (x_left+x_right)/2.0;
			if(abs(x_right-x_left) > OMEGA_STEP) poles.push_back(pole);
		}
	}

	return poles;
	// if(poles.size()==roots.size())
	// {
	// 	return poles;
	// }
	// else
	// {
  //   // roots.resize(poles.size(),0.0);
	// 	return poles;
	// }
}

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
	cout << "sigma = " << sigma.col(2).transpose() << endl << endl;
  // for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
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

	string filename = "debug_f_"+ to_string(int(temperature*100))+".dat";
	ofstream debug_f(filename);
	vector < pair<double,double>> fx_list = get_fx_list(temperature, vt, spa_eivals, debug_f);
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "get_fx_list ");

	vector <double> poles = find_poles_old(fx_list);			
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "find_poles ");

	vector <double> roots = find_roots_old(fx_list,temperature);
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	show_time(begin_ms, end_ms, "roots calc ");

	cout << "Poles " << endl;
	for(int it=0; it<poles.size(); it++) cout << poles.at(it) << endl;
  cout << endl;
	cout << "Roots " << endl;
	for(int it=0; it<roots.size(); it++) cout << roots.at(it) << endl;
  cout << endl;
  
	pair <double, double> free_energies= rpa_free_energy(spa_spectrum.second, roots, poles, temperature);
	cout << temperature << " " << free_energies.first/L << " " << free_energies.second/L << endl;

	cout << "Program is finished." << endl << endl;

	return 0;

}	