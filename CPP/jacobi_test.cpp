#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#include "jacobi.h"
#include "decomposition.h"
#include "operator_matrix.h"

using namespace std;

int main(){
  int myrank, n_procs;
	MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
	
  int Nx = 10,Ny =10,N = Nx*Ny, maxiter=1000;
	double 	Lx = 1.,  Ly = 1., D =1., eps = 10e-10;
	
	decomposition Dc(myrank, n_procs, 1, Nx, Ny);
	operator_matrix A(Nx, Ny, Lx, Ly, D);	
  double *U,*RHS,*Uit,*RHSit;
	U    = (double*) calloc(N,sizeof(double)); 
  RHS  = (double*) calloc(N,sizeof(double));
	for(int i=0;i<N;++i){
		RHS[i] = 1;
		U[i] = 1.1;
	}
	JacobiMethod J(&A,&Dc,RHS,U);
	J.compute(maxiter,eps);
	J.save();	
	MPI_Finalize();
}