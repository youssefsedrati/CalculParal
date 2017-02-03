#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>
#include "tools.h"
#include "conjgrad.h"
#include "cdt_bords.h"
#include "assert.h"
#include "decomposition.h"
#include "operator_matrix.h"


using namespace std;

//void nloc(int *i, int *j, int n, int Nx);
//void matvec(operator_matrix a, decomposition *D, double *Uold, double *U);
//void RightHandSide(int N, int Nx, int M, double dx, double dy,
  //         double Cx, double Cy,double *RHS);

int main(){
  int myrank, nb_procs;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
  
  int Nx = 10,Ny =10,N = Nx*Ny, maxiter=1000;
  double  Lx = 1.,  Ly = 1., D =1., eps = 10e-7;
  decomposition Dc(myrank, nb_procs, 1, Nx, Ny);
  operator_matrix A(Nx, Ny, Lx, Ly, D);
  
    
  double *U,*Uold,*RHS;
  U    = (double*) calloc(N,sizeof(double)); 
  Uold = (double*) calloc(N,sizeof(double));
  RHS  = (double*) calloc(N,sizeof(double));
  //RightHandSide(N, Nx, Ny, dx, dy, Cx, Cy, RHS);
  for(int i=0;i<N;++i){
    RHS[i] = 1;
    U[i] = 1.1;
  }
  /*int *idx = Dc.get_index_global();
  for(int i=0;i<Nx;++i){
    for(int j=0;j<Ny;++j)
      cout << idx[i+j*(Nx)] <<" ";
    cout << endl;
  }*/
  CGmethod C(A, &Dc, RHS, U); // ok
  C.compute(maxiter,eps); // ok
  //C.save(); // pas encore
  
}
/*
int cdt_choisie= 1;
cdt_aux_bords cdt[] = 
{
  {f2, func_zero, func_zero},
  {f1, f1, f1}
};
*/





