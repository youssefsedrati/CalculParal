#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "jacobi.h"
#include "decomposition.h"
#include "operator_matrix.h"

using namespace std;

void fill_RHS_force(decomposition *D,operator_matrix *A, double *RHS);
void fill_RHS_NeumannBC(decomposition *D,operator_matrix *A, double *RHS);

int main(){
  int myRank, n_procs;
	MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
	
  int Nx = 20,Ny =20,N = Nx*Ny, maxiter=10000;
	double 	Lx = 1.,  Ly = 1., D =1., eps = 1e-10, t1,t2;
	bool NeumannBC = false;
	
	if(!myRank) t1=MPI_Wtime(); 
	decomposition Dc(myRank, n_procs, 1, Nx, Ny);
	operator_matrix A(Nx, Ny, Lx, Ly, D,NeumannBC);	
  double *U,*RHS;
	U    = (double*) calloc(N,sizeof(double)); 
  RHS  = (double*) calloc(N,sizeof(double));
	for(int i=0;i<N;++i){
		U[i] = 1.1;
	}
	fill_RHS_force(&Dc,&A,RHS);
	//fill_RHS_NeumannBC(&Dc,&A,RHS);
	JacobiMethod J(&A,&Dc,RHS,U);
	J.compute(maxiter,eps,NeumannBC);
	J.save();	
	MPI_Barrier(MPI_COMM_WORLD);
	if(!myRank) {
		t2=MPI_Wtime(); 
		cout << "time elapsed: " << t2-t1 << endl;
	}
	MPI_Finalize();
}

void fill_RHS_force(decomposition *D,operator_matrix *A, double *RHS){
	int Nx=D->get_Nx(), Ny=D->get_Ny();
	for(int i=0;i<Ny;++i)
		for(int j=0;j<Nx;++j){
			double x = (double)j/(Nx-1), y = (double)i/(Ny-1);
			RHS[j+Nx*i] //= 2* ( x*(1-x) + y*(1-y) );
				//= sin(x)+cos(y);
				= 1;
		}
}

void fill_RHS_NeumannBC(decomposition *D,operator_matrix *A, double *RHS){
	int Nx=D->get_Nx(), Ny=D->get_Ny();
	double dx=A->dx(), dy=A->dy();
	for(int i=1;i<Nx-1;++i){
		double x = (double)i/(Nx-1), y = 0;
		RHS[i] -= //2*(sin(x) + cos(y))/dx;
			0;//2/dx;
		y = 1;
		RHS[i+Nx*(Nx-1)]-= //2*(sin(x) + cos(y))/dx;
			0;//2/dx;
	}
	for(int i=1;i<Ny-1;++i){
		double x = 0, y =(double)i/(Ny-1);
		RHS[i*Nx] -= //2*(sin(x) + cos(y))/dy;
			2/dy;
		x = 1;
		RHS[(i+1)*Nx-1]-= //2*(sin(x) + cos(y))/dy;
			2/dy;
	}
}