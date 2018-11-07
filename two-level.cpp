#include <iostream>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include <complex>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double U;
double t=1;

bool zheev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  double* w = new  double [N];
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = int(A.size())*4;
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [3*LDA];

  zheev_( &eigenvec_choice, &Uchar, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

bool zgeev_cpp(MatrixXcd& A, vector<double>& lambda, char eigenvec_choice='N')
{  
  int N = A.cols();
  int LDA = A.outerStride();
  int INFO = 0;
  __complex__ double* w = new __complex__ double [N];
  __complex__ double* vl;
  __complex__ double* vr;
  char Nchar = 'N';
  char Vchar = 'V';
  char Uchar = 'U';
  int LWORK = pow(2, N);
  __complex__ double* WORK= new __complex__ double [LWORK];
  double* RWORK = new double [LWORK];
  
  zgeev_(&Nchar, &eigenvec_choice, &N, reinterpret_cast <__complex__ double*> (A.data()), &LDA, w, vl, &LDA, vr, &LDA, WORK, &LWORK, RWORK, &INFO );

  lambda.clear();
  for(int i=0; i<N; i++) lambda.push_back(__real__ w[i]);

  delete[] w; delete[] RWORK; delete[] WORK;
  return INFO==0;
}

vector <double> stdEigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda; 
  if(diagonalization_routine(A,lambda,'N')) return lambda;
}

VectorXd Eigenvalues(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'N'))
 	{
		Map<ArrayXd> b(lambda.data(),lambda.size());
  	return b;
	}
}

MatrixXcd Eigenvectors(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
	std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return A; 
}

pair<MatrixXcd, vector<double>> stdEigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V')) return make_pair(A,lambda);
}

pair<MatrixXcd, VectorXd> Eigenspectrum(MatrixXcd A, bool (*diagonalization_routine)(MatrixXcd&, vector <double>&, char)=&zheev_cpp)
{
  std::vector<double> lambda;
  if(diagonalization_routine(A,lambda,'V'))
 	{
    Map<ArrayXd> b(lambda.data(),lambda.size());
    return make_pair(A,b);
	}
}

inline double fermi_fn(double e_minus_mu, double T) {return (isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}


int main()
{
   MatrixXd H  = MatrixXd::Zero(4,4);
   MatrixXd V1 = MatrixXd::Zero(4,4);
   MatrixXd V2 = MatrixXd::Zero(4,4);

   Vector2d m; m << 1 , 1;

   U = 5;

   H << -U*m(0)/2, -t, 0, 0,
         -t, -U*m(1)/2, 0, 0,
         0, 0, U*m(0)/2, -t,
         0, 0, -t, U*m(1)/2;
   
   cout << H << endl << endl;

  pair<MatrixXcd, VectorXd> spa_spectrum = Eigenspectrum(H);
  MatrixXd u = spa_spectrum.first.real();
  VectorXd hf = spa_spectrum.second;
  cout << spa_spectrum.second.transpose() << endl << endl << spa_spectrum.first << endl;

  V1 << -U/2*m(0), 0, 0, 0,
         0, 0 , 0, 0,
         0, 0, U/2*m(0), 0,
         0, 0, 0, 0;

  V2 <<  0, 0, 0, 0,
         0, -U/2*m(1) , 0, 0,
         0, 0, 0, 0,
         0, 0, 0, U/2*m(0);

  MatrixXd v1t = u.inverse()*V1*u;
  MatrixXd v2t = u.inverse()*V2*u;

  cout << v1t << endl << endl << v2t << endl;
  double expr = 1.0;
  double T= 0.01;
  double Omega;

  for(int i=0;  i<4; i++)
  {
     for(int j=0; j<4; j++)
     {
        expr += U*v1t(i,j)*v1t(j,i)*fermi_fn(hf(i),T)-fermi_fn(hf(j),T)/(hf(i)-hf(j)+ Omega);
     }
  }
}