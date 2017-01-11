#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>

#include "jacobi.h"
#include "cdt_bords.h"
#include "assert.h"
#include "decomposition.h"
#include "operator_matrix.h"
#include "comm_ctrl.h"

using namespace std;

void nloc(int *i, int *j, int n, int Nx);

void RightHandSide(int N, int Nx, int M, double dx, double dy,
				   double Cx, double Cy,double *RHS);

int main(){
  int myrank, nb_procs;
	MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
	
  int Nx = 10,Ny =10,N = Nx*Ny, maxiter=1000;
	double 	Lx = 1.,  Ly = 1., D =1., eps = 10e-7;
	decomposition Dc(myrank, nb_procs, 1, Nx, Ny);
	operator_matrix A(Nx, Ny, Lx, Ly, D);
	comm_ctrl C(decomposition *d, double *rhs, double *rhs_up, double *u, double *U_up, 
						int bottomrank, int toprank, int leftrank, int rightrank);
					
  double *U,*RHS;
	U    = (double*) calloc(N,sizeof(double)); 
  RHS  = (double*) calloc(N,sizeof(double));
	for(int i=0;i<N;++i){
		RHS[i] = 1;
		U[i] = 1.1;
	}
	JacobiMethod J(&A,&Dc,RHS,U);
	J.compute(maxiter,eps);
	J.save();
	
}

int cdt_choisie= 1;
cdt_aux_bords cdt[] = 
{
  {f2, func_zero, func_zero},
  {f1, f1, f1}
};

void RightHandSide(int N, int Nx, int M, double dx, double dy, double Cx, double Cy,double *RHS)
{
  int i,j,l,k;
  double posx,posy;

  //M = N/Nx ; /* # de lignes */
  assert(N/Nx == M);
  double (*f)(double, double, double) = cdt[cdt_choisie].f;
  double (*g)(double, double, double) = cdt[cdt_choisie].g;
  double (*h)(double, double, double) = cdt[cdt_choisie].h;

  for( i = 1; i<=N; i++ ){
    nloc(&j,&k,i,Nx);
    posx = k*dx;
    posy = j*dy;
    RHS[i] = f(posx,posy,0.0);
  }

  /* premiere ligne condition de bord du bas */
  
  for( i = 1; i<= Nx; i++ ){
  nloc(&j,&k,i,Nx);
  posx = k*dx;
  posy = j*dy;
  RHS[i] = RHS[i]-g(posx,0.0,0.0)*Cy;
      }

  /* derniere ligne condition de bord du haut */
  l = 1;
  for( i = N-Nx+1;i<=N;i++ ){ 
      nloc(&j,&k,i,Nx);
      posx = k*dx;
      posy = j*dy;
      RHS[i] =RHS[i]-g(posx,1.0,0.0)*Cy;
    }

  /* Bords droit et gauche */
    /*Ligne du bas*/
  RHS[1]  = RHS[1]  -h(0.0,dy,0.0)*Cx;
  RHS[Nx] = RHS[Nx] -h(1.0,dy,0.0)*Cx;

  /*Ligne du milieux*/
  j = 1+Nx;
  for( i = 2; i<= M-1; i++ ){
    nloc(&k,&l,j,Nx);
    RHS[j] = RHS[j] -h(0.0,k*dy,0.0)*Cx;
    RHS[j+Nx-1] = RHS[j+Nx-1] -h(1.0,k*dy,0.0)*Cx;
    j = 1 + (i)*Nx;
  }
  /*ligne du haut*/
  nloc(&k,&l,N,Nx);
  RHS[N-Nx+1] = RHS[N-Nx+1] -h(0.0,k*dy,0.0)*Cx;
  RHS[N] = RHS[N] -h(1.0,k*dy,0.0)*Cx;
}

void nloc(int *i, int *j, int n, int Nx)
{
  int q,r;
  q = n/Nx;
  r = n - q*Nx;  
  if ( r == 0 ){
    *i = q;
    *j = Nx;
  }
  else{
    *i = 1+q;
    *j = r;
  }
  return;
}