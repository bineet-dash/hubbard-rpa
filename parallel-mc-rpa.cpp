#include "alhassid_parallel.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1;
double U_prime=2;
int L=4;


using namespace std::chrono;
using std::cout;
typedef pair <double,double> pdd;

void greens_sigma_generate(MatrixXd& suggested_sigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_sigma(lattice_index,2) *= -1;
}

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) no of sweeps.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  int no_sweeps = atoi(argv[3]);
  int N_therm = 0.5*no_sweeps;
  int N_meas = no_sweeps-N_therm;

  int initial_exp = -3;
  int final_exp = 0;
  double final_temp = 10*pow(10,final_exp);
  milliseconds begin_ms, end_ms;
  long idum = time(NULL);

  MatrixXd sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  for(int i=0; i<L; i++)  greens_sigma_generate(sigma, i, idum);
  MatrixXd suggested_sigma = sigma;
  MatrixXcd H0 = construct_h0();
  MatrixXcd Id = MatrixXcd::Identity(2*L,2*L);

  MPI_Init(NULL, NULL);

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma) + U_prime*L/4*Id;
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);

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

	double pspa_f_parallel = get_parallel_pspa_F(u, spa_eivals, temperature);

  string filename, latticedata;
  latticedata = "_U="+to_string(int(U_prime))+"_size="+to_string(L)+"_sweeps="+to_string(no_sweeps);
  // filename="data/spin_arrangement"+current_time_str()+latticedata+".nb"; ofstream outfile_spinarr(filename);
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  filename="data/m_length_tda_"+ current_time_str()+latticedata+".txt"; ofstream outfile_mlength(filename);
  filename="data/rpa_results_"+current_time_str()+latticedata+".txt"; ofstream outfile_freeenergy(filename);
  // filename="data/mcdetails"+current_time_str()+latticedata+".txt"; ofstream outfile_mcdetails(filename);
  cout << "==============================\n"<< "filename is: " << filename << "\n========================\n";

  // for(int j=final_exp; j>=initial_exp; j--)
  // {
  //   for(double i=10; i>=2; i-=1)
  //   {
  //     double temperature = i*pow(10,j);
  double decrement = 0.05;
  for(double temperature = 2.0; temperature >= 0.01; temperature -= decrement)
  {
    decrement = (temperature > 0.9)?  0.05:0.01;
    for(int sweep=0; sweep<N_therm; sweep++)
    {
      for(int lattice_index=0; lattice_index<L; lattice_index++)
      {
        greens_sigma_generate(suggested_sigma,lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma)+ U_prime*L/4*Id;
        pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
        pdd suggested_free_energies = get_spa_pspa_F(suggested_spa_spectrum.first.real(), suggested_spa_spectrum.second, temperature); 

        double move_prob = exp(-(suggested_free_energies.second-free_energies.second)/temperature);
        double uniform_rv = ran0(&idum);

        if(uniform_rv <= move_prob)
        {
          sigma = suggested_sigma;
          free_energies = suggested_free_energies;
        }
        else
        {
          suggested_sigma=sigma;
        }
      }
      cout << "\r sweep = " << sweep << " done."; cout.flush();
    }

    double final_free_energy_rpa = 0.0;
    double final_free_energy_spa = 0.0;
    double S_pi = 0.0;

    for(int sweep= N_therm; sweep<no_sweeps; sweep++)
    {
      for(int lattice_index=0; lattice_index<L; lattice_index++)
      {
        greens_sigma_generate(suggested_sigma,lattice_index, idum);
        MatrixXcd suggested_Hspa = H0-U_prime/2*matrixelement_sigmaz(suggested_sigma)+ U_prime*L/4*Id;
        pair<MatrixXcd,VectorXd> suggested_spa_spectrum = Eigenspectrum(suggested_Hspa);
        pdd suggested_free_energies = get_spa_pspa_F(suggested_spa_spectrum.first.real(), suggested_spa_spectrum.second, temperature); 

        double move_prob = exp(-(suggested_free_energies.second-free_energies.second)/temperature);
        double uniform_rv = ran0(&idum);

        if(uniform_rv <= move_prob)
        {
          sigma = suggested_sigma;
          free_energies = suggested_free_energies;
        }
        else
        {
          suggested_sigma=sigma;
        }
      }

      final_free_energy_spa += free_energies.first;
      final_free_energy_rpa += free_energies.second; 
      S_pi += get_spi(sigma);
      cout << "\r sweep = " << sweep << " done."; cout.flush();
    }

    outfile_mlength << temperature <<  " " << sigma.col(2).transpose() << endl;
    outfile_freeenergy << temperature << " " << final_free_energy_spa/double(N_meas) << " " << final_free_energy_rpa/double(N_meas) 
                        << " " << S_pi/double(N_meas) << endl;

    cout << "\rtemperature = " << temperature << " done."; cout.flush();
  }

  cout << endl;
  end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  // show_time(begin_ms, end_ms,"MC calculation");
  // spinarrangement_Mathematica_output(sigma,outfile_spinarr);
  // outfile_spinarr.close();

  // outfile_mcdetails.close();
  outfile_mlength.close();
  outfile_freeenergy.close();
  return 0;
}
