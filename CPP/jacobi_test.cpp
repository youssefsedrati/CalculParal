#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "RHS.h"
#include "jacobi.h"
#include "decomposition.h"
#include "operator_matrix.h"

using namespace std;

int main(){
  int myRank, n_procs;
	MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
	
  int Nx = 20,Ny =20,N = Nx*Ny, maxiter=50000, overlap = 0;
	double 	Lx = 1.,  Ly = 1., D =1., eps = 1e-9*N, t1,t2;
	bool NeumannBC = false;
	
	if(!myRank) t1=MPI_Wtime(); 
	decomposition Dc(myRank, n_procs, 1, Nx, Ny,overlap);
	operator_matrix A(Nx, Ny, Lx, Ly, D,NeumannBC);	
  double *U,*RHS;
	U    = (double*) calloc(N,sizeof(double)); 
  RHS  = (double*) calloc(N,sizeof(double));
	fill_RHS_force(&Dc,&A,U,&f1);
	fill_RHS_force(&Dc,&A,RHS,&one);
	//fill_RHS_NeumannBC(&Dc,&A,RHS);
	fill_RHS_DirichletBC(&Dc,&A,RHS,&null);
	JacobiMethod J(&A,&Dc,RHS,U);
	J.compute(maxiter,eps,NeumannBC);
	//J.save();
	J.save_gnuplot();	
	MPI_Barrier(MPI_COMM_WORLD);
	if(!myRank) {
		t2=MPI_Wtime(); 
		cout << "time elapsed: " << t2-t1 << endl;
	}
	MPI_Finalize();
}
