#include <iostream>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include <complex>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double epsilon = 1;
int g = 2;
double chi;

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


MatrixXd construct_H(double sigma0)
{
	MatrixXd H = MatrixXd::Zero(2*g,2*g);
	for(int i=0; i<2*g; i++) H(i,i) = pow(-1,i)*epsilon;
	for(int i=0; i<2*g-1; i+=2) H(i,i+1) = -chi*sigma0;
	for(int i=1; i<2*g; i+=2) H(i,i-1) = -chi*sigma0;
	return H;
}

inline double del(int a1, int a2){return (a1==a2)?1:0;}

double comm_expectation(int i, int j, int v1, int v2, int w1, int w2, int k, int l, MatrixXd U)
{
	double result = 0.0;
	for(int x1=0; x1<2*g; x1++)
	{
		for(int x2=0; x2<2*g; x2++)
		{
			for(int x3=0; x3<2*g; x3++)
			{
				for(int x4=0; x4<2*g; x4++)
				{
					result += U(x1,v1)*U(x2,v2)*U(x3,w1)*U(x4,w2)*(del(i,x1)*del(j,x3)-del(i,x3)*del(j,x1))*(del(l,x2)*del(k,x4)-del(k,x2)*del(l,x4)) ;
				}					
			}			
		}		
	}
	cout << result << endl;
	return result;
}

double lipkin_ijkl(int i, int j, int k, int l, MatrixXd U)
{
	double result = 0.0;
	for(int v=0; v<g; v++)
	{
		for(int w=0; w<g; w++)
		{
			result += comm_expectation(i,j,2*v,2*v+1,2*w,2*w+1,k,l,U) + comm_expectation(i,j,2*v,2*v+1,2*w+1,2*w,k,l,U) +
								comm_expectation(i,j,2*v+1,2*v,2*w,2*w+1,k,l,U) + comm_expectation(i,j,2*v+1,2*v,2*w+1,2*w,k,l,U);
		}
	}
	return -0.5*chi*result;
}

inline double fermi_fn(double e_minus_mu, double T) {return (isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}

MatrixXd construct_A(VectorXd hf, MatrixXd U, double T)
{
	MatrixXd A = MatrixXd::Zero(g*g,g*g);
	for(int j=0; j<2*g; j+=2) 
	{
		for(int i=1; i<2*g; i+=2)
		{
			int it1 = j/2*g+(i-1)/2;
			for(int l=0; l<2*g; l+=2) 
			{
				for(int k=1; k<2*g; k+=2)
				{
					int it2 = l/2*g+(k-1)/2;
					A(it1,it2) = (hf(i)-hf(j))*del(i,k)*del(j,l)+ lipkin_ijkl(i,l,j,k,U)*(fermi_fn(hf(k),T)-fermi_fn(hf(l),T));
				}
			}
		}
	}
	return A;
}

MatrixXd construct_B(VectorXd hf, MatrixXd U, double T)
{
	MatrixXd B = MatrixXd::Zero(g*g,g*g);
	for(int j=0; j<2*g; j+=2) 
	{
		for(int i=1; i<2*g; i+=2)
		{
			int it1 = j/2*g+(i-1)/2;

			for(int l=0; l<2*g; l+=2) 
			{
				for(int k=1; k<2*g; k+=2)
				{
					int it2 = l/2*g+(k-1)/2;
					B(it1,it2) = lipkin_ijkl(i,k,j,l,U)*(fermi_fn(hf(k),T)-fermi_fn(hf(l),T));
				}
			}
		}
	}
	return B;
}

double squareroot(double x){return sqrt(x);}

int main(int argc, char* argv[])
{
	if(argc!=4) {cerr << "Enter (1) chi (2) beta (3) sigma0\n"; exit(1);}
	chi = atof(argv[1]);
	double beta = atof(argv[2]);
	double sigma0 = atof(argv[3]);

	auto H = construct_H(sigma0);
	VectorXd hf_eivals = Eigenvalues(H);
	MatrixXd U = Eigenvectors(H).real();

	MatrixXd A = construct_A(hf_eivals, U, 1/beta);
	MatrixXd B = construct_B(hf_eivals, U, 1/beta);

	cout << H << endl << endl << hf_eivals.transpose() << endl << endl;
	cout << U << endl << endl;
	cout << A << endl << endl << B << endl << endl;

	MatrixXd D = (A+B)*(A-B);
	cout << Eigenvalues(D).transpose() << endl;
}