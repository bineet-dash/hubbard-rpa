#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "rpa.hpp"
#include <clocale>

using namespace std;
using namespace Eigen;

double t=1;
double U_prime=2;
int L=4;

const char  *ptr = NULL;
const wchar_t uparrow[] = L"\u2191";
const wchar_t downarrow[] = L"\u2193";
const wchar_t rangle[] = L"\u3009";

void vis(int u_tilde)
{
  setlocale(LC_ALL, "");
  if(u_tilde==0)      std::wcout << "1" << uparrow;
  else if(u_tilde==1) std::wcout << "2" << uparrow;
  else if(u_tilde==2) std::wcout << "1" << downarrow;
  else if(u_tilde==3) std::wcout << "2" << downarrow;
}

void vis_basis(int c0, int c1)
{
  freopen(ptr, "w", stdout);
  vis(c0);vis(c1); wcout << " ";
  freopen(ptr, "w", stdout);
}

MatrixXd kprod(MatrixXd m1, MatrixXd m2)
{
  MatrixXd m3(m1.rows()*m2.rows(), m1.cols()*m2.cols());

  for (int i = 0; i < m1.cols(); i++)
  {
    for (int j = 0; j < m1.rows(); j++) 
    {
      m3.block(i*m2.rows(), j*m2.cols(), m2.rows(), m2.cols()) =  m1(i,j)*m2;
    }
  }
  return m3;
}

cd debug_rpa_matrix_elem_A(int y1, int y2, int y3, int y4, double T, VectorXd hf, MatrixXcd V_tilde, bool ch = false)
{
  cd u_matrix_elem = 0.25*V_tilde(2*L*y1+y4, 2*L*y2+y3);  //minj -> v_mjin so 1234 -> v_1423
  cd elem = del(y1,y3)*del(y2,y4)*(hf(y1)-hf(y2))+ u_matrix_elem*(fermi_fn(hf(y3),T)-fermi_fn(hf(y4),T));
  // if(ch) cout << u_matrix_elem << " " << u_matrix_elem*(fermi_fn(hf(y3),T)-fermi_fn(hf(y4),T)) << endl;
  return elem;
}

cd debug_rpa_matrix_elem_B(int y1, int y2, int y3, int y4, double T, VectorXd hf, MatrixXcd V_tilde)
{
  cd u_matrix_elem = V_tilde(2*L*y1+y3, 2*L*y2+y4); //minj -> v_mnij so 1234 -> v_1324
  cd elem = u_matrix_elem*(fermi_fn(hf(y3),T)-fermi_fn(hf(y4),T));
  return elem;
}

VectorXd debug_RPA_eivals(VectorXd spa_eivals_minus_mu, MatrixXcd V_tilde, double T)
{
  int N = spa_eivals_minus_mu.size();
  MatrixXcd A = MatrixXcd::Zero(N*N/4,N*N/4);
  MatrixXcd B = MatrixXcd::Zero(N*N/4,N*N/4);

  for(int j=0; j<N/2; j++)
  {
    for(int i=N/2; i<N; i++)
    {
      for(int l=0; l<N/2; l++)
      {
        for(int k=N/2; k<N; k++)
        {
          int it1 = (i-L)*L+j;
          int it2 = (k-L)*L+l; // A_matrix 
          A(it1,it2) = debug_rpa_matrix_elem_A(i, j, k, l, T, spa_eivals_minus_mu, V_tilde);
          B(it1,it2) = debug_rpa_matrix_elem_B(i, j, k, l, T, spa_eivals_minus_mu, V_tilde);
        }
      }
    }
  }

  // cout << "*******************************\n";
  // cout << debug_rpa_matrix_elem_A(3,1,3,1, T,spa_eivals_minus_mu, V_tilde, true) << "\n****************************************\n";

  cout << "A=\n" << A.real() << endl << endl << "B=\n" <<  B.real() << endl << endl;

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

cd debug_U_elem(MatrixXcd u_tilde, int x1, int x2, int x3, int x4);

pair<int,int> get_basis(int x) { int site = x%L; int spin = (x<L)?1:-1; return make_pair(site,spin); }

int main(int argc, char* argv[])
{
  if(argc!=4) {cerr << "Enter (1) lattice size, (2) u_tilde, (3) temp.\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);

  MatrixXd V = MatrixXd::Zero(2*L*2*L,2*L*2*L);

  for(int s1=0; s1<2*L; s1++)
  {
    for(int s2=0; s2<2*L; s2++)
    {
      for(int s3=0; s3<2*L; s3++)
      {
        for(int s4=0; s4<2*L; s4++)
        {
          int it1 = s1*2*L + s2;
          int it2 = s3*2*L + s4;

          pair<int,int> ket1 = get_basis(s1);
          pair<int,int> ket2 = get_basis(s2);
          pair<int,int> ket3 = get_basis(s3);
          pair<int,int> ket4 = get_basis(s4);

          bool same_site = (ket1.first==ket2.first && ket2.first==ket3.first && ket3.first==ket4.first);
          bool up_down_pair = (ket1.second+ket2.second==0) && (ket3.second+ket4.second==0);
          bool minus_U = (ket1.second==ket3.second && ket2.second==ket4.second);

          V(it1, it2) = (same_site && up_down_pair)? 2*(minus_U? -U_prime:U_prime):0.0;
        }
      }
    }
  }
  
  cout.precision(3);


  MatrixXd sigma = MatrixXd::Zero(L,3);
  sigma.col(2) = VectorXd::Constant(L,1);
  for(int i=0; i<L; i++) sigma(i,2) = pow(-1,i);
  MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  MatrixXd u = spa_spectrum.first.real();
  cout << "SPA eivals = " << spa_spectrum.second.transpose() << endl;
  cout << "u = \n" << u << endl << endl;

  MatrixXd u_tilde = MatrixXd::Zero(2*L*2*L,2*L*2*L);
  for(int it1=0;  it1< 4*L*L; it1++) 
  {
    for(int it2=0;  it2< 4*L*L; it2++) 
    {
      int x1 = it1/(2*L); int x2 = it1%(2*L);
      int x3 = it2/(2*L); int x4 = it2%(2*L);
     
      u_tilde(it1,it2) = u(x1,x3)*u(x2,x4);
    }
  }

  cout << "u_tilde= " << endl <<  u_tilde << endl << endl;

  MatrixXd Id = MatrixXd::Identity(H_spa.cols(), H_spa.rows());
  MatrixXd H_spa2 = kprod(H_spa.real(), Id)+kprod(Id, H_spa.real());

  cout << "H_tilde diagonalized=\n" <<  (u_tilde.adjoint()*H_spa2*u_tilde).unaryExpr(&filter_d) << endl << endl;

  cout << "V (in position basis) = \n";
  cout << "    ";
  for(int it=0;  it< 4*L*L; it++) 
  {
    vis_basis(it/(2*L),it%(2*L));
    cout << " ";
  }
  cout << endl;

  for(int it1=0;  it1< 4*L*L; it1++) 
  {
    vis_basis(it1/(2*L),it1%(2*L)); cout << " ";
    for(int it2=0;  it2< 4*L*L; it2++) 
    {
      // if(V(it1,it2)!=0) cout << V(it1,it2) << "  ";
      // else cout << "   ";
      double next_elem = (it2 < V.cols()-1)? V(it1,it2+1):0.0;

      if(next_elem ==0)
        cout << V(it1,it2) << "    ";
      else if(next_elem < 0)
        cout << V(it1,it2) << "  ";
      else 
        cout << V(it1,it2) << "   ";
    }
    cout << endl;
  }
  cout << endl;

  MatrixXd V_tilde = (u_tilde.adjoint()*V*u_tilde).unaryExpr(&filter_d);
  cout << "V_tilde = " << endl;
  for(int it=0;  it< 4*L*L; it++) 
  {
    cout << "   " << it/(2*L) << it%(2*L) << " ";
  }
  cout << endl << endl << V_tilde << endl << endl;

  MatrixXd check_v = MatrixXd::Zero(2*L*2*L,2*L*2*L);
  for(int x1=0; x1<2*L; x1++)
  {
    for(int x2=0; x2<2*L; x2++)
    {
      for(int x3=0; x3<2*L; x3++)
      {
        for(int x4=0; x4<2*L; x4++)
        {
          int it1 = x1*2*L + x2;
          int it2 = x3*2*L + x4;
          check_v(it1,it2) = debug_U_elem(u, x1, x2, x3, x4).real();
        }
      }
    }
  }
  cout << "Check_V= \n" << check_v << endl << endl;
  cout << "Check_V in position basis= \n" << u_tilde*check_v*u_tilde.adjoint() << endl << endl;
  cout << "RPA eivals from V_tilde = " << debug_RPA_eivals(spa_spectrum.second, V_tilde, temperature).transpose() << endl << endl;
  cout << "RPA eivals from check_V = " << debug_RPA_eivals(spa_spectrum.second, check_v, temperature).transpose() << endl << endl;
}

cd debug_U_elem(MatrixXcd u_tilde, int x1, int x2, int x3, int x4)
{
  cd res = 0;
  for(int s=0; s<L; s++) 
  {
    res += conj(u_tilde(x1,s)*u_tilde(x2,s+L))*u_tilde(x3,s)*u_tilde(x4,s+L);
  }
  return (U_prime)*res;
}

  // MatrixXd H_tilde = MatrixXd::Zero(2*L*2*L,2*L*2*L);
  // H_tilde.block(0,0,2*L,2*L) = H_spa.real();
  // H_tilde.block(2*L,2*L,2*L,2*L) = H_spa.real();


