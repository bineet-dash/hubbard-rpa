#include "alhassid_rpa.hpp"

double t=1;
double U_prime=2;
int L=4;

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U, (3) config_int.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
	int config_int = atoi(argv[3]);

	int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);

	MatrixXd sigma = MatrixXd::Zero(L,3);
	string filename = "comparison/spa_vs_pspa_L"+to_string(L)+"_U"+to_string(int(U_prime))+"_cfg"+to_string(config_int)+".dat";
	ofstream dataout(filename);

	sigma.col(2) = 	get_field(config_int);// VectorXd::Constant(L,1); 
	cout << "config: " << sigma.col(2).transpose() <<  "\n"; //cout.flush();

	MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma) + U_prime*L/4*MatrixXcd::Identity(2*L,2*L);
	pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
	MatrixXd u = spa_spectrum.first.real();
	VectorXd spa_eivals = spa_spectrum.second;

	for(int j=final_exp; j>=initial_exp; j--)
  {
    for(double i=10; i>=2; i-=1)
    {
			double temperature = i*pow(10,j);
			pair <double, double> free_energies= get_free_energy_shortcut(u,spa_eivals,temperature);
			dataout << temperature << " " << free_energies.first << " " << free_energies.second << endl;
			std::cout << temperature << " done!" << "\r"; cout.flush();
		}
	}

	std::cout << "Program is finished." << endl << endl;

	return 0;

}	
