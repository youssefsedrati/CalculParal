#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


void Jacobi(int maxiter,int eps,double Aii,double Cx,double Cy,int Nx,int N,double *RHS,double *U, double *Uold)
{
	int myrank, nb_procs;
  	int l, i, M;
  	double residu, drl, dwl, alpha, beta;
  	double *r, *kappa, *d, *W;
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
 	MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
 	MPI_Status s1, s2;
  	MPI_Request r1, r2, r3, r4;
  	double *vect;
  	//initialisation des variables pour le calcul.
  	int i;
  	for(i = 0 ; i<fin ; i++)
  	{	
  		vect[i] = U[i];
  	}

  	// Blocks until all processes have reached this routine.
  	MPI_Barrier(MPI_COMM_WORLD);

}