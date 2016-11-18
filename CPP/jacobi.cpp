#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "tools.h"
#include "jacobi.h"

// constructor & destructor
JacobiMethod::JacobiMethod(double aii,double cx,double cy,int nx,int n,double *rhs,double *u){
	Aii = aii; Cx = cx; Cy = cy; Nx = nx; N = n; RHS = rhs; U = u; Ny = N/Nx;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
}

JacobiMethod::~JacobiMethod(){
	cleanup();
}
// public
void JacobiMethod::compute(int itermax, double e){
	init(itermax,e);
	compute_iterate();
	cleanup();	
}
// private
void JacobiMethod::init(int itermax, double e){
	iterMax = itermax; eps = e;
  iter = 0;
  dist = 1;
  Uit = (double*) calloc(N+1,sizeof(double)); 
  d   = (double*) calloc(N+1,sizeof(double));
  W   = (double*) calloc(N+1,sizeof(double));
	init_MPI();
	init_sys();
}

void JacobiMethod::init_MPI(){
  firstline = myrank*Ny/nb_procs+1;
  lastline = (myrank+1)*Ny/nb_procs;
  if(myrank == nb_procs-1) lastline = Ny;

  first = (firstline-1) * Nx + 1;
  last = lastline * Nx;	
}

void JacobiMethod::init_sys(){
	for(int i = 0 ; i<=N ; i++) 
		Uit[i] = U[i];	
}

void JacobiMethod::compute_iterate(){
	while( (iter<iterMax)&&(dist>eps*eps) ){
		if(iter%2) compute_alternate_update(Uit,U);
		else compute_alternate_update(U,Uit);
		dist= dist_squared(N,U,Uit);
		iter+=1;
	}
}

void JacobiMethod::compute_alternate_update(double* Uold, double* Unew){
	for(int i=0;i<=N;++i){
		sigma=0;
		if(i!=0) sigma=+Cx*Uold[i];
		if(i!=N) sigma=+Cx*Uold[i];
		if(i>=Ny) sigma=+Cy*Uold[i];
		if(i<=N-Ny) sigma=+Cy*Uold[i];
		Unew[i]= (RHS[i]-sigma)/Aii;
	}
}

void JacobiMethod::compute_gen_sol(){
	if(iter%2)
		for(int i=0;i<=N;++i) U[i]= Uit[i];
		MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0){
    printf("Jacobi method terminated after"
           " %d iterations, squared_dist= %0.12f\n", iter-1,dist);
  }
}

void JacobiMethod::cleanup(){
	if(Uit) free(Uit);
	if(d) free(d);
	if(W) free(W);
}