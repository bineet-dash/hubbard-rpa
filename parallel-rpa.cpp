#include "alhassid_rpa.hpp"
#include "alhassid_parallel.hpp"

double t=1;
double U_prime=2;
int L=4;

milliseconds begin_ms, end_ms;


int main(int argc, char* argv[])
{
	MPI_Init(NULL, NULL);

  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) temp, (4) config_int .\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);
	int config_int = atoi(argv[4]);


	MatrixXd sigma = MatrixXd::Zero(L,3);
	sigma.col(2) = get_field(config_int);
  MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  MatrixXd u = spa_spectrum.first.real();
  VectorXd spa_eivals = spa_spectrum.second;
	
	// cout << spa_eivals.transpose() << endl;

	
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	if(process_rank==0)
	{
		begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		double pspa_f = get_pspa_F(u,spa_eivals, temperature);
		end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		cout << pspa_f << endl;
		cout << (end_ms-begin_ms).count() << " miliseconds" << endl << endl;
	}

	begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	double pspa_f_parallel = get_pspa_F(u,spa_eivals, temperature);
	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());

	if(process_rank==0)
	{
		cout << "pspa_F " <<  pspa_f_parallel << endl;
		cout << (end_ms-begin_ms).count() << " miliseconds" << endl << endl;
	}	
		

	
	MPI_Finalize();
	return 0;
}	

/* 	ofstream no_root("no-root.dat");
	int matsubara_r = 0;
	double final_det_r = 0.0;

	int r_max = int(abs( (spa_eivals(spa_eivals.size()-1)-spa_eivals(0))/temperature )) ; //omega_max = (2r_max+1)*pi*T= \delta_ij_max
	cout << "r_max = " << r_max << endl;
	do
	{
		double omega_r = (2* matsubara_r +1)*M_PI*temperature;
		double det_r = log(rpa_det_cd(omega_r, vt, spa_eivals, temperature));
		no_root << matsubara_r << " " << det_r << endl;
		final_det_r -= det_r;
		matsubara_r ++; 
	}
	while(matsubara_r < r_max); 


	cout << -temperature*final_det_r << endl;
	cout << (end_ms-begin_ms).count() << endl; */

	// pair <double, double> free_energies= rpa_free_energy(spa_spectrum.second, roots, poles, temperature);
	// cout << temperature << " " << free_energies.first/L << " " << free_energies.second/L << endl;