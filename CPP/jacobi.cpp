#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "jacobi.h"
#include "decomposition.h"
#include "operator_matrix.h"
#include "comm_ctrl.h"

/* implementation of the iterative Jacobi Method.
	 the solution vectors are updated alternatingly;
	 each process uses its decomposition (ie. subdomain) to operate on;
	 additive schwarz is realized by transmitting Dirichlet data via comm_ctrl;
*/

using namespace std;

// constructor & destructor
JacobiMethod::JacobiMethod(operator_matrix *a,decomposition *dc,
				double *rhs, double *u){
	A = a; D = dc; RHS = rhs; U = u;
	myN  = D->get_myN();
	myNx = D->get_myNx();
	myNy = D->get_myNy();
	N  = D->get_N();
	Nx = D->get_Nx();
	Ny = D->get_Ny();
  Uit  = (double*) calloc(N,sizeof(double)); // solution is alternatingly stored in U and Uit
	RHSit= (double*) malloc(N*sizeof(double)); // RHS updated via comm_ctrl
	C = new comm_ctrl(D,A,RHS,RHSit);
}

JacobiMethod::~JacobiMethod(){
	//cleanup();
	delete C;
}
// public
void JacobiMethod::compute(int itermax, double e, bool neumannBC){
	init(itermax,e, neumannBC);
	compute_iterate();	
	compute_gen_sol();
}

void JacobiMethod::save(){
	//if(myRank!=0) return;
	cout << "#" << myRank << ". saving.\n";
	std::string filename = "Jacobi_test_";
	filename.append(to_string(myRank)); filename.append(".data");
	std::ofstream ofs(filename,std::ofstream::out);
	if(!ofs) return;
	ofs << Nx << " "<< Ny <<endl;
	for(int i=Ny-1;i>=0;--i){
		for(int j=0;j<Nx;++j)
			ofs << U[j+Nx*i] << " ";
		ofs << endl;
	}
	ofs.close();
}

void JacobiMethod::save_gnuplot(){ 
  if(myRank!=0) return; 
  cout << "#" << myRank << ". saving.\n"; 
  std::string filename = "Jacobi_gnuplot.dat"; 
  std::ofstream ofs("Jacobi_gnuplot.dat",std::ofstream::out); 
  if(!ofs) return; 
  for(int i=0;i<Ny;++i){ 
    for(int j=0;j<Nx;++j){ 
      double x= (double)(j+1)/(Nx+1), y= (double)(i+1)/(Ny+1); 
      ofs << x << " " << y << " " << U[j+Nx*i] << endl; 
    } 
    ofs << endl; 
  } 
  ofs.close(); 
} 

// private
void JacobiMethod::init(int itermax, double e, bool neumannBC){
	NeumannBC = neumannBC;
	iterMax = itermax; eps = e;
  iter = dist_squared = 1;
	init_MPI();
	init_sys();
}

void JacobiMethod::init_MPI(){	
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nOfProcs);
}

void JacobiMethod::init_sys(){
	for(int i=0;i<N;++i){
		Uit[i] = 0;	
		RHSit[i] = RHS[i];
	}
}

void JacobiMethod::compute_iterate(){	
	double eps_sq = eps*eps;
	while( (iter<=iterMax)&&(dist_squared>eps_sq) ){
		C->receive(); 
		compute_iterate_local();
		C->send(U);
		C->cumulate_dist_squared(&dist_squared);
	}
}

void JacobiMethod::compute_iterate_local(){
	int it=0;
	double eps_ = eps*eps;
	if(D->get_N_procs()>1) eps_*= D->get_N_procs()*100;
	while( (it<iterMax)&&(dist_squared>eps_) ){
		if(iter%2){
			compute_alternate_update(Uit,U);
		}else{
			compute_alternate_update(U,Uit);
		}
		compute_dist_squared();
		it++;
		iter+=1;
	}
}

void JacobiMethod::compute_alternate_update(double* U, double* Unew){
	compute_alternate_update_inner(U, Unew);
	compute_alternate_update_bottom(U, Unew);
	compute_alternate_update_left(U, Unew);
	compute_alternate_update_right(U, Unew);
	compute_alternate_update_top(U, Unew);
}

void JacobiMethod::compute_alternate_update_inner(double* U, double* Unew){
	int *idx = D->get_index_global_inner();
	for(int j=0;j<D->get_myNinner();++j){
		int i = idx[j];
		sigma = A->Cx()*( U[i-1]+U[i+1] )
					 +A->Cy()*( U[i-Nx]+U[i+Nx] );
		Unew[i]= (RHSit[i]-sigma)/A->Aii();
	}	
}

void JacobiMethod::compute_alternate_update_left(double* U, double* Unew){
	int i,*idx = D->get_index_global_left();
	if(NeumannBC && D->get_myRank_x()==0){
		i = idx[0];
		Unew[i] = (RHSit[i]-A->Cx()*U[i+1]-A->Cy()*U[i+Nx])/A->Aii()*2;
		i = idx[myNy-1];	
		Unew[i] = (RHSit[i]-A->Cx()*U[i+1]-A->Cy()*U[i-Nx])/A->Aii()*2;
	}else{
		i = idx[0];
		Unew[i] = (RHSit[i]-A->Cx()*U[i+1]-A->Cy()*U[i+Nx])/A->Aii();
		i = idx[myNy-1];	
		Unew[i] = (RHSit[i]-A->Cx()*U[i+1]-A->Cy()*U[i-Nx])/A->Aii();
	}
	for(int j=1;j<myNy-1;++j){
		i = idx[j];
		if(NeumannBC && D->get_myRank_x()==0)
			sigma = 2*A->Cx()* U[i+1] +A->Cy()*( U[i-Nx]+U[i+Nx] );
		else
			sigma = A->Cx()* U[i+1] +A->Cy()*( U[i-Nx]+U[i+Nx] );
		Unew[i]= (RHSit[i]-sigma)/A->Aii();
	}
}

void JacobiMethod::compute_alternate_update_right(double* U, double* Unew){
	int i,*idx = D->get_index_global_right();
	if(NeumannBC && D->get_myRank_x()==D->get_N_procs_x()-1){
		i = idx[0];
		Unew[i] = (RHSit[i]-A->Cx()*U[i-1]-A->Cy()*U[i+Nx])/A->Aii()*2;
		i = idx[myNy-1];	
		Unew[i] = (RHSit[i]-A->Cx()*U[i-1]-A->Cy()*U[i-Nx])/A->Aii()*2;
	}else{
		i = idx[0];
		Unew[i] = (RHSit[i]-A->Cx()*U[i-1]-A->Cy()*U[i+Nx])/A->Aii();
		i = idx[myNy-1];	
		Unew[i] = (RHSit[i]-A->Cx()*U[i-1]-A->Cy()*U[i-Nx])/A->Aii();
	}
	for(int j=1;j<myNy-1;++j){
		i = idx[j];
		if(NeumannBC && D->get_myRank_x()==D->get_N_procs_x()-1)
			sigma = 2*A->Cx()* U[i-1] +A->Cy()*( U[i-Nx]+U[i+Nx] );
		else
			sigma = A->Cx()* U[i-1] +A->Cy()*( U[i-Nx]+U[i+Nx] );
		Unew[i]= (RHSit[i]-sigma)/A->Aii();
	}	
}

void JacobiMethod::compute_alternate_update_bottom(double* U, double* Unew){
	int *idx = D->get_index_global_bottom();
	for(int j=1;j<myNx-1;++j){
		int i = idx[j];
		if(NeumannBC && D->get_myRank_y()==0)
			sigma = A->Cx()*( U[i-1]+U[i+1] ) +2*A->Cy()* U[i+Nx];
		else
			sigma = A->Cx()*( U[i-1]+U[i+1] ) +A->Cy()* U[i+Nx];
		Unew[i]= (RHSit[i]-sigma)/A->Aii();
	}
}

void JacobiMethod::compute_alternate_update_top(double* U, double* Unew){
	int *idx = D->get_index_global_top();
	for(int j=1;j<myNx-1;++j){
		int i = idx[j];		
		if(NeumannBC && D->get_myRank_y()==D->get_N_procs_y()-1)
			sigma = A->Cx()*( U[i-1]+U[i+1] ) +2*A->Cy()* U[i-Nx];
		else
			sigma = A->Cx()*( U[i-1]+U[i+1] ) +A->Cy()* U[i-Nx];
		Unew[i]= (RHSit[i]-sigma)/A->Aii();
	}	
}

void JacobiMethod::compute_dist_squared(){
	double Di; dist_squared=0;
	int *idx = D->get_index_global();
	for(int i=0;i<myN;++i){
		int j = idx[i];
		Di= U[j]-Uit[j];
		dist_squared+= Di*Di;
	}
}

void JacobiMethod::compute_gen_sol(){
	for(int i=0;i<N;++i) U[i]= Uit[i];
	MPI_Barrier(MPI_COMM_WORLD);
	C->compile_solution(U);		
  if(myRank == 0){
    printf("Jacobi method terminated after"
           " %d iterations, dist= %0.31f\n", iter-1,sqrt(dist_squared));
  }
}

void JacobiMethod::cleanup(){
	if(Uit) free(Uit);
}
