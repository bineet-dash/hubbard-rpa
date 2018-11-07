#include "rpa.hpp"

double t=1;
double U_prime=2;
int L=4;
MatrixXcd U;

VectorXd debug_RPA_eivals(VectorXd spa_eivals_minus_mu, MatrixXcd U, double T)
{
  int N = spa_eivals_minus_mu.size();
  MatrixXcd A = MatrixXcd::Zero(N*N/4,N*N/4);
  MatrixXcd B = MatrixXcd::Zero(N*N/4,N*N/4);

  for(int j=0; j<N/2; j++)
  {
    for(int i=N/2; i<N; i++)
    {
      int it1 = j*L+(i-L);
      for(int l=0; l<N/2; l++)
      {
        for(int k=N/2; k<N; k++)
        {
          int it2 = l*L+(k-L); // A_matrix 
          A(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu, U);
        }
      }
      for(int l=N/2; l<N; l++)
      {
        for(int k=0; k<N/2; k++)
        {
          int it2 = k+(l-L)*L;
          B(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu, U);
        }
      }
    }
  }

  cout << A << endl << endl << B << endl << endl;

  MatrixXcd D = (A+B)*(A-B);
  VectorXd raw_rpa_eivals = Eigenvalues(D).unaryExpr(&filter_d);

  if(raw_rpa_eivals.minCoeff() >=0 )
  {
    VectorXd rpa_eivals = raw_rpa_eivals.unaryExpr(&squareroot);
    return rpa_eivals;
  }
  else
  {
    cerr << "Imaginary RPA Eigenvalues! D-matrix eigenvalues are: " << endl << raw_rpa_eivals.transpose() << endl;
    return VectorXd::Zero(raw_rpa_eivals.size());
  }
}


int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) U and (3) temperature.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);

  MatrixXd sigma = MatrixXd::Zero(L,3);
  // sigma.col(2) = VectorXd::Constant(L,1);
  for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
  MatrixXcd H0 = construct_h0();

  MatrixXcd H_spa = H0 - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);

  cout << H_spa.real() << endl << endl << spa_spectrum.second.transpose() << endl << endl << spa_spectrum.first.real() << endl << endl;
  cout << sigma.col(2).transpose() << endl;
  cout << debug_RPA_eivals(spa_spectrum.second, spa_spectrum.first, temperature).transpose() << endl;



}