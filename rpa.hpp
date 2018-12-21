#ifndef _RPA_HPP_INCLUDED
#define _RPA_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <lapacke.h>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>

using namespace std;
using namespace Eigen;

typedef std::complex <double> cd;

extern double t;
extern double U_prime;
extern int L;

inline double del(int a1, int a2){return (a1==a2)?1:0;}
inline cd jn(cd z){return conj(z);}
inline pair<int,int> mi(int index){return make_pair(int(index/L), index%L+L);}
inline double Sqr(double x){return x*x;}
inline double gs_energy(VectorXd hf_eivals) {return hf_eivals.block(0,0,hf_eivals.size()/2,1).sum();}
inline double gs_energy(vector<double> v) {return accumulate(v.begin(), v.begin()+v.size()/2, 0.00);}
inline cd filter_cd(cd x){return (abs(x)<1e-4)?0.0:x;}
inline double filter_d(double x) {return (abs(x)<1e-6)?0.0:x;}
inline double filter_tol_d(double x, double tolerance=1e-4) {return (abs(x)<tolerance)?0.0:x;}
inline double fermi_fn(double e_minus_mu, double T) {return (isinf(exp(e_minus_mu/T)))? 0: 1/(exp(e_minus_mu/T)+1);}


#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
double ran0(long *idum)
{
   long  k;
   double ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

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

MatrixXcd construct_h0(void)
{
  MatrixXcd Mc = MatrixXcd::Zero(2*L,2*L);
  for(int row=0;row<2*L-1; row++) Mc(row,row+1)=Mc(row+1,row)=-t;
  Mc(L-1,0)=Mc(0,L-1)=-t; //PBC
  Mc(2*L-1,L)=Mc(L,2*L-1)= -t; //PBC
  Mc(L,L-1)= Mc(L-1,L)=0;
  return Mc;
}

MatrixXcd matrixelement_sigmax(MatrixXd randsigma)
{
  MatrixXcd Mcx = MatrixXcd::Zero(2*L,2*L);
  for(int row=0; row<L; row++)
  {
    Mcx(row,row+L) = cd(randsigma(row,0),0);
    Mcx(row+L,row) = cd(randsigma(row,0),0);
  }
  return Mcx;
}

MatrixXcd matrixelement_sigmay(MatrixXd randsigma)
{
	MatrixXcd Mcy = MatrixXcd::Zero(2*L,2*L);
	for(int row=0; row<L; row++)
	{
		Mcy(row,row+L) = cd(0,-randsigma(row,1));
		Mcy(row+L,row) = cd(0,randsigma(row,1));
	}
	return Mcy;
}

MatrixXcd matrixelement_sigmaz(MatrixXd randsigma)
{
	MatrixXcd Mcz = MatrixXcd::Zero(2*L,2*L);
	for(int row=0; row<L; row++)
		Mcz(row,row)= cd(randsigma(row,2),0);
	for(int row=L; row<2*L; row++)
		Mcz(row, row)=cd(-randsigma(row-L,2),0);
	return Mcz;
}

double get_spi(MatrixXd sigma )
{
  double sq = 0.0;
  for(int i=0; i<L; i++)
  {
    for(int j=0; j<L; j++)
    {
      sq += sigma(i,2)*sigma(j,2)*pow(-1,i-j)/pow(L,2);
    }
  }
  return sq;
}

double get_mu(double temperature, std::vector<double> v)
{
  sort (v.begin(), v.end());
  double bisection_up_lim = v.back();
  double bisection_low_lim = v.front();

  double mu, no_of_electrons; int count=0;
  double epsilon = 0.000001;

  for(; ;)
  {
    no_of_electrons=0;
    mu = 0.5*(bisection_low_lim+bisection_up_lim);

    for(auto it = v.begin(); it!= v.end(); it++)
    {
      double fermi_func = 1/(exp((*it-mu)/temperature)+1);
      no_of_electrons += fermi_func;
    }
    if(abs(no_of_electrons-L) < epsilon)
    {
      return mu; break;
    }
    else if(no_of_electrons > L+epsilon)
    {
       if(abs(bisection_up_lim-v.front())<0.001){return mu; break;}
       else {bisection_up_lim=mu;}
    }
    else if(no_of_electrons < L-epsilon)
    {bisection_low_lim=mu;}
  }
}

double get_mu(double temperature, VectorXd v)
{
  vector<double> stdv (v.data(),v.data()+v.size());
  return get_mu(temperature, stdv);
}

double spa_free_energy(VectorXd spa_eivals, double temperature)
{
  double free_energy = 0; double ekt =0;
  double mu = get_mu(temperature, spa_eivals);

  for(auto it=0; it!= spa_eivals.size(); it++)
  {
    ekt = (spa_eivals[it]-mu)/temperature;
    if(!isinf(exp(-ekt))) free_energy += -temperature*log(1+exp(-ekt));
    else  free_energy += (spa_eivals[it]-mu);
  }
  return free_energy/L+mu;
}


double spa_internal_energy(MatrixXcd Mc, double temperature)
{
  std::vector<double> eigenvalues;
  zheev_cpp(Mc, eigenvalues, 'N');
  sort(eigenvalues.begin(),eigenvalues.end());
  double mu = get_mu(temperature, eigenvalues);

  double internal_energy=0.0 ; double e_min = eigenvalues.front();
  for(auto it=eigenvalues.begin(); it!= eigenvalues.end(); it++)
  {
    internal_energy += (*it)/(exp((*it-mu)/temperature)+1);
  }
  return internal_energy;
}

cd U_elem(MatrixXcd U, int x1, int x2, int x3, int x4)
{
  cd res = 0;
  for(int s=0; s<L; s++) 
  {
    res += conj(U(s,x1)*U(s+L,x2))*U(s,x3)*U(s+L,x4);
    // res += conj(U(s,l)*U(s+L,i)-U(s,i)*U(s+L,l))*(U(s,k)*U(s+L,j)-U(s,j)*U(s+L,k));
  }
  return (U_prime)*res;
}

cd rpa_matrix_elem(int i, int j, int k, int l, double T, VectorXd hf, MatrixXcd U)
{
  cd u_iljk = 2.0*U_elem(U.transpose(),i,l,j,k);
  cout << u_iljk.real() << endl;
  cd elem = del(i,k)*del(j,l)*(hf(i)-hf(j)) + u_iljk*(fermi_fn(hf(k),T)-fermi_fn(hf(l),T));
  return elem;
}

MatrixXcd construct_RPA(VectorXd spa_eivals, MatrixXcd U, double T)
{
  double mu = get_mu(T, spa_eivals);
  VectorXd spa_eivals_minus_mu = spa_eivals.array()-mu;
  int N = spa_eivals.size();
  
  MatrixXcd H_rpa = MatrixXcd::Zero(N*N/2,N*N/2);

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
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu, U);
          H_rpa(it1+2*L, it2+2*L) = -conj(H_rpa(it1,it2));
        }
      }

      for(int l=N/2; l<N; l++)
      {
        for(int k=0; k<N/2; k++)
        {
          int it2 = k+(l-L)*L+ L*L;
          H_rpa(it1,it2) = rpa_matrix_elem(i, j, k, l, T, spa_eivals_minus_mu,U);
          H_rpa(it1+2*L, it2-2*L) = -conj(H_rpa(it1,it2));
        }
      }
    }
  }
  return H_rpa;
}

inline double squareroot(double x){return sqrt(x);}

VectorXd get_RPA_eivals(VectorXd spa_eivals_minus_mu, MatrixXcd U, double T)
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

inline double log_sinh(double x) { return abs(x)+log(1-exp(-2*abs(x))) ; } //1/2 factor not needed

double rpa_free_energy(VectorXd spa_eivals, MatrixXcd U, double T)
{
  double mu = get_mu(T, spa_eivals);
  spa_eivals = spa_eivals.array()-mu;
  VectorXd rpa_eivals = get_RPA_eivals(spa_eivals,U,T);

  double beta = 1/T;
  double spa_part = 0.0, rpa_part_nom = 0.0, rpa_part_denom = 0.0;

  for(int i=0; i<spa_eivals.size(); i++) 
    spa_part += (-beta*spa_eivals(i)>4.0)?  -beta*spa_eivals(i):log(1+exp(-beta*spa_eivals(i)));

  for(int i=0; i<rpa_eivals.size(); i++)
  {
    rpa_part_denom += (rpa_eivals(i)==0)? log(beta/2): (log_sinh(beta*rpa_eivals(i)/2) - log(abs(rpa_eivals(i))) );
  }

  for(int j=0; j<spa_eivals.size(); j++)
  {
    for(int i=0; i<j; i++)
    {
      double Deltaij = filter_d(spa_eivals(i)-spa_eivals(j));
      if(Deltaij==0) continue;
      else
      {
        rpa_part_nom += ( log_sinh(beta*Deltaij/2) - log(abs(Deltaij)) );
                // cout << Deltaij << ", " << rpa_part_nom << " -> ";
      }
    }
  }
  // cout << endl;

  // cout.precision(7); //cout << "eivals = " << rpa_eivals.transpose() << " \n\n";
  // cout <<  setw(10) << -T*rpa_part_nom/L << " " << setw(10) << T*rpa_part_denom/L << " " << setw(10) << -T*spa_part/L << "  ";

  return -T*(rpa_part_nom-rpa_part_denom+spa_part)/L + mu;
}

string current_time_str(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];
  time (&rawtime); 
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%S-%M-%I-%Y-%m-%d",timeinfo);
  string str(buffer);
  return str;
}

void spinarrangement_Mathematica_output(MatrixXd M, ofstream& fout)
{
  double lattice_separation = 1.0;
  fout << "Show[Graphics3D[{" << endl;
  for(int i=0; i< M.rows(); i++)
  {
    fout << "Arrow[{{" << lattice_separation*i << ", 0, 0}, {" << M(i,0)+lattice_separation*i << ","  << M(i,1) << "," << M(i,2) << "}}]";
    if(i!=M.rows()-1) fout << ",\n";
  }
  fout <<"}] ]" << endl << endl;
}

#endif


/* double get_moveProbability(VectorXd spa, VectorXd spa_new, VectorXd rpa, VectorXd rpa_new, double T)
{
  double beta = 1/T;
  VectorXd delta_spa = spa_new - spa;
  VectorXd delta_rpa = rpa_new - rpa;

  double rpa_part = 0.0;
  for(int i=0; i<spa.size(); i++)
  {
    for(int j=0; j<spa.size(); j++)
    {
      if(i==j) continue;
      if(spa(i)==spa(j))
        rpa_part += -log(spa_new(i)-spa_new(j))+ beta*(delta_spa(i)-delta_spa(j)) + log(1-exp(-2*beta*(spa_new(i)-spa_new(j))));

      else if(spa_new(i)==spa_new(j)) 
        rpa_part += log(spa(i)-spa(j))+ beta*(delta_spa(i)-delta_spa(j)) - log(1-exp(-2*beta*(spa(i)-spa(j))));

      else
        rpa_part += log(spa(i)-spa(j))-log(spa_new(i)-spa_new(j))+ beta*(delta_spa(i)-delta_spa(j))
                + log(1-exp(-2*beta*(spa_new(i)-spa_new(j)))) - log(1-exp(-2*beta*(spa(i)-spa(j))));
    }
  }
  for(int i=0; i<rpa.size(); i++)
    rpa_part += log(rpa_new(i)/rpa(i)) - beta*delta_rpa(i) + log(1-exp(-2*beta*(rpa(i)))) - log(1-exp(-2*beta*rpa_new(i)));
  
  double spa_part = 0.0;
  for(int i=0; i<rpa.size(); i++)
    spa_part += (-beta*spa_new(i) > 4.0)?-beta*spa_new(i):log(1+exp(-beta*spa_new(i)) ) - (-beta*spa(i) > 4.0)?-beta*spa(i):log(1+exp(-beta*spa(i))) ;
  
  return exp(rpa_part+spa_part);
}
 */