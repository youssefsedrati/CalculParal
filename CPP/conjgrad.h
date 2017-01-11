#ifndef CONJGRAD_H
#define CONJGRAD_H

#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"


class CGmethod{
public:
	CGmethod();
	CGmethod(operator_matrix a, decomposition *dc, double *RHS,double *U);
	~CGmethod();
	void compute(int iterMax, double eps);
private:
// system variables
	operator_matrix A;
	decomposition *D;
	double  eps,residue, drl, dwl, alpha, beta;
	double *RHS = NULL, *U = NULL, *r = NULL, *kappa = NULL, *d = NULL, *W = NULL;
	int iterMax,iter, Nx,Ny,N;
// MPI variables
  int myrank, nb_procs,fstline,lstline,fst,lst,deb,fin;

  MPI_Request r1,r2,r3,r4;
  MPI_Status s1,s2;
// functions
	void init(int itermax, double e);
	void init_sys();
	void init_MPI();
	void compute_iterate();
	void compute_comm_main();
	void compute_calc_pre();
	void compute_calc_main();
	void compute_gen_sol();
	void cleanup();
};


#endif
