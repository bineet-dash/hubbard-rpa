#include "rpa.hpp"
#include <cstring>
#include <chrono>
#include <cstdlib>

double t=1;
double U_prime=2;
int L=4;
MatrixXcd U;

using namespace std::chrono;

void greens_sigma_generate(MatrixXi& suggested_sigma, int lattice_index, long & idum)
{
  if(ran0(&idum)<=0.5) suggested_sigma(lattice_index,2) *= -1;
}

VectorXd inttobin(int theValue)
{
  VectorXd v(L);
  for (int i = 0; i < L; ++i)  v(i) = theValue & (1 << i) ? 1 : 0;
  return v;
}

VectorXd get_field(int i)
{
  VectorXd raw = inttobin(i);
  for(int i=0; i<raw.size(); i++) raw(i) = (raw(i)==0)?-1:1;
  return raw;
}

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) temperature.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);

  MatrixXd sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  MatrixXcd H0 = construct_h0();
  
  int i_min = 0, i_max = 0;
  for(int it=0; it < L; it++) i_max += pow(2,it);

  {
    MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma);
    pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
    double free_energy = rpa_free_energy(spa_spectrum.second, spa_spectrum.first, temperature);
    cout << "Ferro: " << free_energy/L << "  " << " " << endl;
  }
 
  {
    for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
    MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma);
    pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
    double free_energy = rpa_free_energy(spa_spectrum.second, spa_spectrum.first, temperature);
    cout << "Antiferro: " << free_energy/L << "  " << " " << endl;
  }
  // MatrixXcd H_rpa = construct_RPA(spa_spectrum.second, spa_spectrum.first, temperature);
  // cout.precision(3);

  // // for(int i=0; i<N*N/4; i++) cout << "(" << 
  // cout << endl << H_rpa.unaryExpr(&filter).real() << endl << endl;

  // MatrixXcd A = H_rpa.block(0, 0, H_rpa.rows()/2, H_rpa.cols()/2); 
  // MatrixXcd B = H_rpa.block(0, H_rpa.cols()/2, H_rpa.rows()/2, H_rpa.cols()/2);
  // MatrixXcd D = (A+B)*(A-B);

  // cout << Eigenvalues(D).unaryExpr(&squareroot).transpose() << endl << endl;
  // cout << Eigenvalues(H_rpa, &zgeev_cpp).transpose() << endl << endl;
}