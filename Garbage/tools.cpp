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

void matvec(operator_matrix a, decomposition *D, double *Uold, double *U){
  int     i,j,k;
  int myrank;
  int nb_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int fst_line = myrank* D->get_myNy()/nb_procs+1;
  int lst_line = (myrank+1)* D->get_myNy()/nb_procs;
  if(myrank == nb_procs-1)
    lst_line = D->get_myNy();

  int fst_elt = (fst_line-1) * D->get_myNx() + 1;
  i = fst_elt;
  // Je suppose pas de problème avant cette ligne
  if(fst_line == 1){
/*Premier bloc*/
    U[1] = a.Aii()*Uold[1] + a.Cx()*Uold[2] + a.Cy()*Uold[1+ D->get_myNx()];
    i = 1;
    for( j = 1;j<= D->get_myNx()-2;j++ ){
      i = 1+j;
      U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i-1] + a.Cx()*Uold[i+1] + a.Cy()*Uold[i+D->get_myNx()];
    }
    U[1+(D->get_myNx()-1)] = a.Aii()*Uold[1+(D->get_myNx()-1)] + a.Cx()*Uold[1+(D->get_myNx()-1)-1] + a.Cy()*Uold[1+(D->get_myNx()-1)+D->get_myNx()];
    i = 1 + (D->get_myNx()-1);
  }
/*bloc general, il y a m-2 blocs generaux */
  for( k = fst_line; k<= MIN(D->get_myNy()-1, lst_line); k++){
    if(k!=1){ /* First line already done */
/*Premiere ligne*/
      i = (k-1)* D->get_myNx()+1;
      U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i+1] + a.Cy()*Uold[i- D->get_myNx()] + a.Cy()*Uold[i+ D->get_myNx()];
/*ligne generale*/
      for( j = 1;j<= D->get_myNx() -2;j++){
        i = i+1;
        U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i-1] + a.Cx()*Uold[i+1] + a.Cy()*Uold[i+ D->get_myNx()] + a.Cy()*Uold[i- D->get_myNx()];
      }
/*Derniere ligne*/
      i = i + 1;
      U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i-1] + a.Cy()*Uold[i- D->get_myNx()] + a.Cy()*Uold[i+ D->get_myNx()];
    }
  }
  i = i+1;
  if(lst_line == D->get_myNy()){
/*Dernier bloc*/
    U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i+1] + a.Cy()*Uold[i- D->get_myNx()];
    for( j = 1;j<= D->get_myNx()-2;j++){
      i = i+1;
      U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i+1] + a.Cx()*Uold[i-1] + a.Cy()*Uold[i- D->get_myNx()];
    }
    i = i+1;
    U[i] = a.Aii()*Uold[i] + a.Cx()*Uold[i-1] + a.Cy()*Uold[i- D->get_myNx()];
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
