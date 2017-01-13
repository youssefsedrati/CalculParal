#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include "decomposition.h"
#include "comm_ctrl.h"

using namespace std;

int main(){
	int myRank, nOfProcs, nOfProcs_x=2, Nx = 10,Ny =10,N = Nx*Ny;
	MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nOfProcs);
	
	decomposition Dc(myRank, nOfProcs, nOfProcs_x, Nx, Ny);
	
	double 	Lx = 1.,  Ly = 1., D =1.;
	operator_matrix A(Nx, Ny, Lx, Ly, D);	

  double *U,*RHS,*Uit,*RHSit;
	U    = (double*) malloc(N*sizeof(double)); 
  RHS  = (double*) malloc(N*sizeof(double)); 
  Uit  = (double*) malloc(N*sizeof(double)); 
	RHSit= (double*) malloc(N*sizeof(double));
	for(int i=0;i<N;++i){
		RHS[i] = 1;
		U[i] = 1.1;
	}

	comm_ctrl C(&Dc,&A,RHS,RHSit,U,Uit);

	
	MPI_Finalize();
	return 0;
}