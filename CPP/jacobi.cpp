#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "tools.h"
#include "jacobi.h"
#include "decomposition.h"
#include "operator_matrix.h"
#include "comm_ctrl.h"

using namespace std;

// constructor & destructor
JacobiMethod::JacobiMethod(operator_matrix *a,decomposition *dc,double *rhs,double *u,comm_ctrl *c){
	A = a; D = dc; RHS = rhs; U = u; C = c;
	N  = D->get_myN();
	Nx = D->get_myNx();
	Ny = D->get_myNy(); 
  Uit = (double*) malloc(N*sizeof(double)); 
	RHSit =(double*) malloc(N*sizeof(double));
}

JacobiMethod::~JacobiMethod(){
	//cleanup();
}
// public
void JacobiMethod::compute(int itermax, double e){
	init(itermax,e);
	compute_iterate();
	compute_gen_sol();
}

void JacobiMethod::save(){
	FILEOUT = ofstream("Jacobi_test_sol.data",ios::out);
	if(!FILEOUT) return;
	FILEOUT << Nx << " "<< Ny<<endl;
	for(int i=0;i<Nx;++i){
		for(int j=0;j<Ny;++j)
			FILEOUT << U[i+Ny*j] << " ";
		FILEOUT << endl;
	}
	FILEOUT.close();
}

// private
void JacobiMethod::init(int itermax, double e){
	iterMax = itermax; eps = e;
  iter = 1;
  dist = 1;
	init_MPI();
	init_sys();
}

void JacobiMethod::init_MPI(){	
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
}

void JacobiMethod::init_sys(){
	for(int i = 0 ; i<=N ; i++){
		Uit[i] = U[i];	
		RHSit[i] = RHS[i];
	}
}

void JacobiMethod::compute_iterate(){	
	while( (iter<iterMax)&&(dist>eps*eps) ){
		C->receive(); // eventually wait until all processors have sent
		if(iter%2) compute_alternate_update(Uit,U);
		else compute_alternate_update(U,Uit);
		C->send(); // eventually wait until all processors have received
		dist= dist_squared(N,U,Uit);
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
		sigma = A->Cx()*U[i-1]+A->Cx()*U[i+1]+A->Cy()*U[i-Ny]+A->Cy()*U[i+Ny];
		Unew[i]= (RHS[i]-sigma)/A->Aii();
	}	
}

void JacobiMethod::compute_alternate_update_left(double* U, double* Unew){
	int i,*idx = D->get_index_global_left();
	i = idx[0];
	Unew[i] = (RHS[i]-A->Cx()*U[i+1]-A->Cy()*U[i+Ny])/A->Aii();
	i = idx[D->get_myNy()-1];	
	Unew[i] = (RHS[i]-A->Cx()*U[i+1]-A->Cy()*U[i-Ny])/A->Aii();
	for(int j=1;j<D->get_myNy()-1;++j){
		i = idx[j];
		sigma = A->Cx()*U[i+1]+A->Cy()*U[i-Ny]+A->Cy()*U[i+Ny];
		Unew[i]= (RHS[i]-sigma)/A->Aii();
	}
}

void JacobiMethod::compute_alternate_update_right(double* U, double* Unew){
	int i,*idx = D->get_index_global_right();
	i = idx[0];
	Unew[i] = (RHS[i]-A->Cx()*U[i-1]-A->Cy()*U[i+Ny])/A->Aii();
	i = idx[D->get_myNy()-1];	
	Unew[i] = (RHS[i]-A->Cx()*U[i-1]-A->Cy()*U[i-Ny])/A->Aii();
	for(int j=1;j<D->get_myNy()-1;++j){
		i = idx[j];
		sigma = A->Cx()*U[i-1]+A->Cy()*U[i-Ny]+A->Cy()*U[i+Ny];
		Unew[i]= (RHS[i]-sigma)/A->Aii();
	}	
}

void JacobiMethod::compute_alternate_update_bottom(double* U, double* Unew){
	int *idx = D->get_index_global_bottom();
	for(int j=1;j<D->get_myNx()-1;++j){
		int i = idx[j];
		sigma = A->Cx()*U[i-1]+A->Cx()*U[i+1]+A->Cy()*U[i+Ny];
		Unew[i]= (RHS[i]-sigma)/A->Aii();
	}
}

void JacobiMethod::compute_alternate_update_top(double* U, double* Unew){
	int *idx = D->get_index_global_top();
	for(int j=1;j<D->get_myNx()-1;++j){
		int i = idx[j];
		sigma = A->Cx()*U[i-1]+A->Cx()*U[i+1]+A->Cy()*U[i-Ny];
		Unew[i]= (RHS[i]-sigma)/A->Aii();
	}	
}

double JacobiMethod::dist_squared(int N,double* U,double* V){
	double D = 0, Di;
	int i;
	for(i=0;i<=N;++i){
		Di= U[i]-V[i];
		D+= Di*Di;
	}
	return D;
}

void JacobiMethod::compute_gen_sol(){
	if(iter%2)
		for(int i=0;i<=N;++i) U[i]= Uit[i];
	MPI_Barrier(MPI_COMM_WORLD);
	C->compile_solution();
  if(myrank == 0){
    printf("Jacobi method terminated after"
           " %d iterations, squared_dist= %0.31f\n", iter-1,dist);
  }
}

void JacobiMethod::cleanup(){
	if(Uit) free(Uit);
}