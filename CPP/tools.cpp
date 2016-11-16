#include "tools.h"
#include "cdt_bords.h"

#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>

int cdt_choisie= 1;
cdt_aux_bords cdt[] = 
{
  {f2, func_zero, func_zero},
  {f1, f1, f1}
};


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

void matvec(double Aii,double Cx,double Cy,int Nx,int Ny,double *Uold,double *U){
  int     i,j,k;
  int myrank;
  int nb_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int fst_line = myrank*Ny/nb_procs+1;
  int lst_line = (myrank+1)*Ny/nb_procs;
  if(myrank == nb_procs-1)
    lst_line = Ny;

  int fst_elt = (fst_line-1) * Nx + 1;
  i = fst_elt;

  if(fst_line == 1){
/*Premier bloc*/
    U[1] = Aii*Uold[1] + Cx*Uold[2] + Cy*Uold[1+Nx];
    i = 1;
    for( j = 1;j<= Nx-2;j++ ){
      i = 1+j;
      U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cx*Uold[i+1] + Cy*Uold[i+Nx];
    }
    U[1+(Nx-1)] = Aii*Uold[1+(Nx-1)] + Cx*Uold[1+(Nx-1)-1] + Cy*Uold[1+(Nx-1)+Nx];
    i = 1 + (Nx-1);
  }
/*bloc general, il y a m-2 blocs generaux */
  for( k = fst_line; k<= MIN(Ny-1, lst_line); k++){
    if(k!=1){ /* First line already done */
/*Premiere ligne*/
      i = (k-1)*Nx+1;
      U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cy*Uold[i-Nx] + Cy*Uold[i+Nx] ;
/*ligne generale*/
      for( j = 1;j<= Nx-2;j++){
        i = i+1;
        U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cx*Uold[i+1] + Cy*Uold[i+Nx] + Cy*Uold[i-Nx];
      }
/*Derniere ligne*/
      i = i + 1;
      U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cy*Uold[i-Nx] + Cy*Uold[i+Nx];
    }
  }
  i = i+1;
  if(lst_line == Ny){
/*Dernier bloc*/
    U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cy*Uold[i-Nx];
    for( j = 1;j<= Nx-2;j++){
      i = i+1;
      U[i] = Aii*Uold[i] + Cx*Uold[i+1] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
    }
    i = i+1;
    U[i] = Aii*Uold[i] + Cx*Uold[i-1] + Cy*Uold[i-Nx];
  }
}

double dist_squared(int N,double* U,double* V){
	double D = 0, Di;
	int i;
	for(i=0;i<=N;++i){
		Di= U[i]-V[i];
		D+= Di*Di;
	}
	return D;
}
