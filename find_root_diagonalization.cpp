#include "alhassid_rpa.hpp"
#include <chrono>

using namespace std::chrono;

void show_time(milliseconds begin_ms, milliseconds end_ms, string s)
{
   long int t = (end_ms.count()-begin_ms.count())/1000;
    if (t<=60)
    { cout <<  s << " took " << t << " seconds." << endl; }
    else if (t>60 && t<3600)
    {
      int minutes=int(t/60); int seconds= t%minutes;
      cout << s << " took " << minutes << " minute and " << seconds << " seconds." << endl;
    }
    else if(t>=3600)
    {
      int hrs= int(t/3600); int minutes=int((t-3600*hrs)/60); int seconds= int(t-3600*hrs-60*minutes);
      cout << s << " took " << hrs << " hour, " << minutes << " minutes and " << seconds << " seconds. ";
    }
    else
    {cout << s << " took " << t << "time. Wrong t received.\n"; }
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

double t=1;
double U_prime=2;
int L=4;

milliseconds begin_ms, end_ms;

double u_iljk(vector <MatrixXd> vt, int i, int l, int j, int k)
{
	double sum=0.0;
	for(int alpha=0; alpha<vt.size(); alpha++)
	{
		sum -= U_prime/2*vt.at(alpha)(i,j)*vt.at(alpha)(l,k);
		// cerr << vt.at(alpha)(i,j) << " " << vt.at(alpha)(l,k) << " " << sum << endl;
	}
	return sum;
}

int main(int argc, char* argv[])
{
  if(argc!=5) {cerr << "Enter (1) lattice size, (2) U, (3) temp, (4) config_int .\n"; exit(1);}
  L = atoi(argv[1]);
  U_prime = atof(argv[2]);
  double temperature = atof(argv[3]);
	int config_int = atoi(argv[4]);

	vector <MatrixXd> v;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		v.push_back(v_i);
	}

  begin_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());


	MatrixXd sigma = MatrixXd::Zero(L,3);
	sigma.col(2) = get_field(config_int);
	cout << "sigma = " << sigma.col(2).transpose() << endl;
  MatrixXcd H_spa = construct_h0() - U_prime/2*matrixelement_sigmaz(sigma);
  pair<MatrixXcd,VectorXd> spa_spectrum = Eigenspectrum(H_spa);
  MatrixXd u = spa_spectrum.first.real();
  VectorXd spa_eivals = spa_spectrum.second;
	cout << "spa_eivals = \n" << spa_eivals.transpose() << endl;

	end_ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
	// show_time(begin_ms, end_ms, "spa calc ");

	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i_transformed = u.adjoint()*(v.at(it))*u;
		vt.push_back(v_i_transformed);
	}

	vector <double> poles = find_poles(vt, spa_eivals, temperature);
	vector <double> roots = find_roots(vt, spa_eivals, temperature, poles);
	cout << endl << "Poles " << endl;
	for(int it=0; it<poles.size(); it++) cout << poles.at(it) << endl;
  cout << endl;
	cout << "Roots " << endl;
	for(int it=0; it<roots.size(); it++) cout << roots.at(it) << endl;
  cout << endl;

	cout.precision(3);

	// for(auto i:vt) cout << i << endl << endl;

	double mu = get_mu(temperature, spa_eivals);
	VectorXd hf_mu = spa_eivals.array()-mu;
	
/* 	cout << "debugging" << endl;
	double debg = u_iljk(vt,0,2,2,0)*(fermi_fn(hf_mu(0),temperature)-fermi_fn(hf_mu(2),temperature));
	cout << "debugging finished" << endl << endl;
	exit(1); */

	MatrixXd rpa = MatrixXd::Zero(4*L*L,4*L*L);
	for(int i=0; i<2*L; i++)
	{
		for(int j=0; j<2*L; j++)
		{
			for(int k=0; k<2*L; k++)
			{
				for(int l=0; l<2*L; l++)
				{
					double it1 = i*2*L+j;
					double it2 = k*2*L+l;
					rpa(it1,it2) = -(hf_mu(i)-hf_mu(j))*del(i,k)*del(j,l) +u_iljk(vt,i,l,j,k)*(fermi_fn(hf_mu(k),temperature)-fermi_fn(hf_mu(l),temperature));
				}
			}
		}	
	}

	for(int i=0; i<2*L; i++)
	{
		double it_remove = 2*L*i+i;
		removeRow(rpa,it_remove);
		removeColumn(rpa, it_remove);
	}

	// cout << "RPA matrix" << endl;
	// for(int k=0; k<2*L; k++)
	// {
	// 	for(int l=0; l<2*L; l++)
	// 	{
	// 		if(k==l) continue;
	// 		cout << "    " << k << l << " ";
	// 	}
	// }
	// cout << endl << rpa.unaryExpr(&filter_d) << endl << endl;
	// cout << Eigenvalues(rpa).unaryExpr(&filter_d).transpose() << endl;

  vector <double> delta;
	// for(int i=0; i<spa_eivals.size(); i++)
	// {
	// 	for(int j=0; j<i; j++)
	// 	{
	// 		double delta_ij = spa_eivals(i)-spa_eivals(j);
	// 		delta.push_back(delta_ij);
	// 	}
	// 	sort(delta.begin(), delta.end());
	// }
	for(int i=0; i<spa_eivals.size()/2; i++)
	{
		for(int j=spa_eivals.size()/2; j< spa_eivals.size(); j++)
		{
			double delta_ij = spa_eivals(i)-spa_eivals(j);
			delta.push_back(delta_ij);
		}
		sort(delta.begin(), delta.end());
	}
	for(int j=0; j<spa_eivals.size()/2; j++)
	{
		for(int i=spa_eivals.size()/2; i < spa_eivals.size(); i++)
		{
			double delta_ij = spa_eivals(i)-spa_eivals(j);
			delta.push_back(delta_ij);
		}
		sort(delta.begin(), delta.end());
	}
	cout << endl << "delta_ij = " << endl;
	for(const auto& i:delta) cout << i << " "; cout << endl;


	MatrixXcd rpa_u = MatrixXcd::Zero(2*L*L, 2*L*L);
	MatrixXcd rpa_del_ij = MatrixXcd::Zero(2*L*L, 2*L*L);

	for(int j=0; j<L; j++)
  {
    for(int i=L; i<2*L; i++)
    {
      int it1 = j*L+(i-L);
      for(int l=0; l<L; l++)
      {
        for(int k=L; k<2*L; k++)
        {
          int it2 = l*L+(k-L); // A_matrix 
          rpa_del_ij(it1,it2) = -(hf_mu(i)-hf_mu(j))*del(i,k)*del(j,l);
          rpa_del_ij(it1+L*L, it2+L*L) = -conj(rpa_del_ij(it1,it2));
        }
      }
    }
  }

  for(int j=0; j<L; j++)
  {
    for(int i=L; i<2*L; i++)
    {
      int it1 = j*L+(i-L);
      for(int l=0; l<L; l++)
      {
        for(int k=L; k<2*L; k++)
        {
          int it2 = l*L+(k-L); // A_matrix 
          rpa_u(it1,it2) = u_iljk(vt,i,l,j,k)*(fermi_fn(hf_mu(k),temperature)-fermi_fn(hf_mu(l),temperature));
          rpa_u(it1+L*L, it2+L*L) = -conj(rpa_u(it1,it2));
        }
      }

      for(int l=L; l<2*L; l++)
      {
        for(int k=0; k<L; k++)
        {
          int it2 = k+(l-L)*L+ L*L;
          rpa_u(it1,it2) = u_iljk(vt,i,l,j,k)*(fermi_fn(hf_mu(k),temperature)-fermi_fn(hf_mu(l),temperature));
          rpa_u(it1+L*L, it2-L*L) = -conj(rpa_u(it1,it2));
        }
      }
    }
  }

	MatrixXcd rpa_old = rpa_u + rpa_del_ij;
	// cout << endl << rpa_u.unaryExpr(&filter_cd).real() << endl << endl;
	cout << "old eivals: " << endl;
	cout << Eigenvalues(rpa_old).block(L*L,0,L*L,1).unaryExpr(&filter_d).transpose() << endl;

	VectorXd raw_eivals = Eigenvalues(rpa_old).transpose();
	MatrixXcd xi_eivec = Eigenvectors(rpa_old);
	// cout << xi_eivec.real().unaryExpr(&filter_d) << endl << endl;

	// cout << "Analysis\n===============\n";

	for(int it=0; it<raw_eivals.size(); it++)
	{
		double Omega = raw_eivals(it);
		// cout << Omega << endl;
		int delet = 1;
		for(int it2=0; it2 < raw_eivals.size(); it2++)
		{
			if(abs(Omega-rpa_old(it2,it2)) > 0.01)
			{
				continue;	
			} 
			else
			{
				// cout << "= del_ij for it2 = " << it2 << ",  i = " << it2%L << ", j = " << it2/L << endl; 
				if(rpa_old.row(it2).dot(xi_eivec.col(it2))== rpa_old(it2,it2)*xi_eivec(it2,it2)) 
				{
					delet*= 1;	
				}
				else
				{
					delet*=0;
				}
				if(delet == 1) 
				{
					raw_eivals(it) = 15690;
				}
			}
		}
	}
	 
	cout << endl << "processed eivals = " << endl << raw_eivals.transpose() << endl << endl;

	// pair <double, double> free_energies= rpa_free_energy(spa_spectrum.second, roots, poles, temperature);
	// cout << temperature << " " << free_energies.first/L << " " << free_energies.second/L << endl;

	return 0;

}	

/* 
{
	auto it=delta .begin();
	while( it!=delta.end()-1)
	{
		if(abs(*(it)-*(it+1)) < 0.01)	delta.erase(it+1);
		else it++;
	}
	for(const auto& i:delta) cout << i << "\n"; cout << endl;
} */