#include "alhassid_rpa.hpp"

double t=1;
double U_prime=2;
int L=4;

int main(int argc, char* argv[])
{
  if(argc!=3) {cerr << "Enter (1) lattice size, (2) U.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);

	int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
	
	MatrixXd sigma = MatrixXd::Zero(L,3);
	string filename = "comparison/total_F_vs_T_L"+ to_string(L)+"_U"+to_string(int(U_prime)) +".dat";
	ofstream dataout(filename);

  int i_min = 0, i_max = 0;
  for(int it=0; it < L; it++) i_max += pow(2,it);

	for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
			double temperature = i*pow(10,j);
			
			double exp_spa_F = 0.0, exp_pspa_F=0.0;
			vector <double> spa_F, pspa_F;

			for(int it= i_min; it<= i_max; it++) 
			{
				sigma.col(2) = get_field(it);
				// cout << "config " << it << " " << sigma.col(2).transpose() <<  "\n"; //cout.flush();
				MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma) + U_prime/4*L*MatrixXcd::Identity(2*L,2*L);
				pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
				MatrixXd u = spa_spectrum.first.real();
				VectorXd spa_eivals = spa_spectrum.second;

				pair <double, double> free_energies= get_free_energy_shortcut(u,spa_eivals,temperature);
				spa_F.push_back(free_energies.first);
				pspa_F.push_back(free_energies.second);
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
			dataout << temperature << " " << spa_energy << " " << pspa_energy << endl;

			std::cout << temperature << " done!" << endl;
		}
	}


	std::cout << "Program is finished." << endl << endl;

	return 0;

}	
