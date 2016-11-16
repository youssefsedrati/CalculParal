#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "tools.h"

void Jacobi(int maxiter,int eps,double Aii,double Cx,double Cy,int Nx,int N,double *RHS,double *U)
{
	int myrank, nb_procs;
	int iter=0, i, j, Ny=N/Nx;
	double dist=1, sigma;
	double *Uit,  *d, *W;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
	MPI_Status s1, s2;
	MPI_Request r1, r2, r3, r4;

	// Init: Communication variables
  Uit = (double*) calloc(N+1,sizeof(double)); 
  d   = (double*) calloc(N+1,sizeof(double));
  W   = (double*) calloc(N+1,sizeof(double));

  int firstline = myrank*Ny/nb_procs+1;
  int lastline = (myrank+1)*Ny/nb_procs;
  if(myrank == nb_procs-1) lastline = Ny;

  int first = (firstline-1) * Nx + 1;
  int last = lastline * Nx;

	// Init: computation variables
	for(i = 0 ; i<=N ; i++) 
		Uit[i] = U[i];
	
	// Computation
	iter = 0;
	while( (iter<maxiter)&&(dist>eps*eps) ){
		// alternate between Uit and U
		if(iter%2){
			for(i=0;i<=N;++i){
				sigma=0;
				if(i!=0) sigma=+Cx*Uit[i];
				if(i!=N) sigma=+Cx*Uit[i];
				if(i>=Ny) sigma=+Cy*Uit[i];
				if(i<=N-Ny) sigma=+Cy*Uit[i];
			}
			U[i]= (RHS[i]-sigma)/Aii;
		}else{
			for(i=0;i<=N;++i){
				sigma=0;
				if(i!=0) sigma=+Cx*U[i];
				if(i!=N) sigma=+Cx*U[i];
				if(i>=Ny) sigma=+Cy*U[i];
				if(i<=N-Ny) sigma=+Cy*U[i];
			}
			Uit[i]= (RHS[i]-sigma)/Aii;
		}
		dist= dist_squared(N,U,Uit);
		iter+=1;
	}
	// make sure U carries the last computation
	if(iter%2)
		for(i=0;i<=N;++i) U[i]= Uit[i];

	// Blocks until all processes have reached this routine.
	MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0){
    printf("Jacobi method terminated after"
           " %d iterations, squared_dist= %0.12f\n", iter-1,dist);
  }
}
