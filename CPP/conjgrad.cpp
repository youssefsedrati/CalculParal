#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "tools.h"
#include "conjgrad.h"
#include "decomposition.h"
#include "operator_matrix.h"

// constructor & destructor
CGmethod::CGmethod(operator_matrix a, decomposition *dc, double *rhs,double *u){
	//Aii = aii; Cx = cx; Cy = cy; Nx = nx; N = n;
	A = a; D = dc; RHS = rhs; U = u; 
	N  = D->get_myN();
	Nx = D->get_myNx();
	Ny = D->get_myNy(); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
}

CGmethod::~CGmethod(){
	//cleanup();
}
// public
void CGmethod::compute(int itermax, double e){
	init(itermax,e);
	compute_iterate();
	cleanup();	
}
// private
void CGmethod::init(int itermax, double e){
	iterMax = itermax; eps = e;
  iter = 0;
  residue = 0.0;
  r     = (double*) calloc(N+1,sizeof(double)); 
  kappa = (double*) calloc(N+1,sizeof(double));
  d     = (double*) calloc(N+1,sizeof(double));
  W     = (double*) calloc(N+1,sizeof(double));
	init_MPI();
	init_sys();
}

void CGmethod::init_MPI(){
  fstline = myrank*Ny/nb_procs+1;
  lstline = (myrank+1)*Ny/nb_procs;
  if(myrank == nb_procs-1)
    lstline = Ny;
  fst = (fstline-1) * Nx + 1;
  lst = lstline * Nx;
  deb = MAX(fst-Nx, 1);
  fin = MIN(lst+Nx,N);	
}

void CGmethod::init_sys(){
	for (int i=deb;i<fin;i++ ){
    kappa[i] = U[i];
  }
  // Problème avec matvec.
  matvec(A,D,kappa,r);

  for(int i=fst;i<=lst;i++){
    r[i]     = r[i] - RHS[i];
    residue = residue + r[i]*r[i];
    d[i]=r[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &residue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void CGmethod::cleanup(){	
	if(r) free(r); 
	if(kappa) free(kappa); 
	if(d) free(d); 
	if(W) free(W);
}

void CGmethod::compute_iterate(){
  while( (iter<=iterMax) && (sqrt(residue) >= eps)){    
		compute_comm_main();
		compute_calc_pre();
		compute_calc_main();
    iter++;
  }
	compute_gen_sol();
}

void CGmethod::compute_comm_main(){
	if(myrank != nb_procs - 1){
		MPI_Isend(&d[lst - Nx + 1], Nx, MPI_DOUBLE,myrank+1, 99, MPI_COMM_WORLD, &r1);
		MPI_Irecv(&d[lst+1], Nx, MPI_DOUBLE,myrank+1, 99, MPI_COMM_WORLD, &r2);
	}
	if(myrank != 0){
		MPI_Isend(&d[fst], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r3);
		MPI_Irecv(&d[fst-Nx], Nx, MPI_DOUBLE, myrank-1, 99, MPI_COMM_WORLD, &r4);
	}
	if(myrank != nb_procs - 1)
	{
		MPI_Wait(&r1, &s1);
		MPI_Wait(&r2, &s1);
	}
	if(myrank != 0)
	{
		MPI_Wait(&r3, &s2);
		MPI_Wait(&r4, &s2);
	}
}

void CGmethod::compute_calc_pre(){

	matvec(A,D,d,W);
	drl = 0.0;
	dwl = 0.0;
	for(int i=fst; i<=lst; i++ ){
		drl += d[i]*r[i];
		dwl += d[i]*W[i];
	}
	MPI_Allreduce(MPI_IN_PLACE, &dwl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &drl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void CGmethod::compute_calc_main(){
	alpha = drl/dwl;
	for(int i=fst; i<=lst; i++ ){
		kappa[i] = kappa[i] - alpha*d[i];
		r[i] = r[i] - alpha*W[i];
	}
	beta = 0.0;
	for(int i=fst;i<=lst;i++){
		beta = beta + (r[i]*r[i]);
	}
	MPI_Allreduce(MPI_IN_PLACE, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	beta = beta / residue;
	residue = 0.0;
	for(int i = fst; i<= lst; i++){
		d[i] = r[i] + beta*d[i];   
		residue    = residue + r[i]*r[i];  
	}
	MPI_Allreduce(MPI_IN_PLACE, &residue, 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}

void CGmethod::compute_gen_sol(){
	for(int i=fst;i<=lst;i++){
    U[i] = kappa[i]; 
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0){
    printf("le Gradient Conjugue a converge en"
           " %d iteration, residue= %0.12f\n", iter,residue);
  }
}