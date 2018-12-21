#include "alhassid_rpa.hpp"

double t=1;
double U_prime=2;
int L=4;

int main(int argc, char* argv[])
{
  if(argc!=3) {cerr << "Enter (1) lattice size, (2) temp.\n"; exit(1);}
  L = atoi(argv[1]);
  double temperature = atof(argv[2]);

	int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
	
	MatrixXd sigma = MatrixXd::Zero(L,3);
	ofstream dataout;
	string filename = "comparison/total_F_vs_U_L"+ to_string(L)+".dat";
	dataout.open(filename);

  int i_min = 0, i_max = 0;
  for(int it=0; it < L; it++) i_max += pow(2,it);


	for(U_prime = 0.0; U_prime <= 25; U_prime += 0.2)
	{
		vector <double> spa_F, pspa_F;

		for(int it= i_min; it<= i_max; it++) 
		{
			sigma.col(2) = get_field(it);
			// cout << "config " << it << " " << sigma.col(2).transpose() <<  "\n"; //cout.flush();
			MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
			pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
			MatrixXd u = spa_spectrum.first.real();
			VectorXd spa_eivals = spa_spectrum.second;

			pair <double, double> free_energies= get_free_energy_shortcut(u,spa_eivals,temperature);
			spa_F.push_back(free_energies.first/L);
			pspa_F.push_back(free_energies.second/L);
		}
		sort(spa_F.begin(),spa_F.end());
		sort(pspa_F.begin(), pspa_F.end());
		// for(auto const& i:spa_F) cout << i << " "; cout << endl;
		// for(auto const& i:pspa_F)  cout << i << " "; cout << endl;

		double remainder_spa_F = 0.0, remainder_pspa_F=0.0;
		for(auto const& i:spa_F) remainder_spa_F += exp(-(i-spa_F[0])/temperature);
		for(auto const& i:pspa_F) remainder_pspa_F += exp(-(i-pspa_F[0])/temperature);
		
		double spa_energy = spa_F[0]-temperature*log(remainder_spa_F);
		double pspa_energy = pspa_F[0]-temperature*log(remainder_pspa_F);
		dataout << U_prime << " " << spa_energy << " " << pspa_energy << endl;

		// dataout << U_prime << " " << -temperature*log(sqrt(U_prime/(4*temperature*M_PI))*exp_spa_F) << " " << -temperature*log(sqrt(U_prime/(4*temperature*M_PI))*exp_pspa_F) << endl;
		// std::cout <<  sqrt(U_prime/(4*temperature*M_PI)) << " " << log( sqrt(U_prime/(4*temperature*M_PI))) << " " << -temperature*log( sqrt(U_prime/(4*temperature*M_PI))) << endl;
		// dataout << temperature << " " << free_energies.first << " " << free_energies.second << endl;
		std::cout << U_prime << " done!" << "\r"; cout.flush();
	}


	std::cout << "Program is finished." << endl << endl;

	return 0;

}	
