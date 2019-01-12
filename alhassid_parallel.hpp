#include "alhassid_rpa.hpp"
#include <mpi.h>

double get_parallel_pspa_F(MatrixXd u, VectorXd hf, double T)
{
	vector <MatrixXd> vt;
	for(int it=0; it<L; it++)
	{
		MatrixXd v_i = MatrixXd::Zero(2*L,2*L);
		v_i(it,it) = 1; v_i(it+L, it+L) = -1;
		MatrixXd v_i_transformed = u.adjoint()*v_i*u;
		vt.push_back(v_i_transformed);
	}

	double mu = get_mu(T, hf);
	VectorXd fermi_hf = VectorXd::Zero(hf.size());
	for(int it=0; it<hf.size(); it++)
	{
		fermi_hf(it) = fermi_fn(hf(it)-mu,T);
	}

	int r_max = int(abs( (hf(hf.size()-1)-hf(0))/T )) ; //omega_max = (2r_max+1)*pi*T= \delta_ij_max
	double final_det_r = 0;

	
	for(int matsubara_r = 0; matsubara_r < 5*r_max; matsubara_r++)
	{
		double omega_r = (2* matsubara_r +1)*M_PI*T;
		RMatrixXcd rpa = MatrixXcd::Zero(L,L);
		int num_procs, rank;

		MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		assert(L*L%num_procs==0);
		int elem_per_proc = L*L/num_procs;
    
    /* if(rank==0)
    {
      cout << "num procs, elem_per procs = " << num_procs << " " << elem_per_proc << endl;
    } */

		cd* localdata = new cd[elem_per_proc];

		MPI_Scatter(rpa.data(),elem_per_proc,MPI_DOUBLE_COMPLEX, localdata, elem_per_proc, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		for(int lin_alpha=0; lin_alpha<elem_per_proc; lin_alpha++)
		{
			int row = (lin_alpha+rank*L)/L; 
			int col = (lin_alpha+rank*L) % L;
      localdata[lin_alpha] = cd(del(row,col),0);
			for(int i=0; i<hf.size(); i++)
			{
				for(int j=0; j<hf.size(); j++)
				{
					cd num = U_prime/2*(vt.at(col))(j,i)*(vt.at(row))(i,j)*(fermi_hf(i)-fermi_hf(j));
					cd denom = cd(hf(i)-hf(j),omega_r);
					localdata[lin_alpha] += num/denom;
				}
			}
		}
		
		MPI_Gather(localdata,elem_per_proc,MPI_DOUBLE_COMPLEX, rpa.data(), elem_per_proc, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		if(rank==0)
		{
			final_det_r += log( real(rpa.determinant()) );
		}
	}
	
	int process_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
	if(process_rank==0)
	{
		return T*final_det_r;
	}
}