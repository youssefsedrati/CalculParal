#ifndef JACOBI_H
#define JACOBI_H

#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"
#include "comm_ctrl.h"

class JacobiMethod{
public:
	JacobiMethod();
	JacobiMethod(operator_matrix *a,decomposition *dc,
				double *rhs, double *u);
	~JacobiMethod();
	void compute(int itermax,double e);
	void save();
private:
// system variables
	operator_matrix *A;
	decomposition *D;
	comm_ctrl *C;
	double eps, dist_squared, sigma;
	double *RHS=NULL,*RHSit=NULL,*U=NULL,*Uit=NULL;
	int iterMax,iter, Nx,Ny,N;
// MPI variables
  int myRank, nOfProcs;
  //MPI_Status s1, s2;
  //MPI_Request r1, r2, r3, r4;
// file operations
	std::ofstream FILEOUT;
// functions
	void init(int itermax, double e);
	void init_sys();
	void init_MPI();
	void compute_iterate();
	void compute_comm_main();
	void compute_alternate_update(double* U1, double* U2);
	void compute_alternate_update_inner(double* U1, double* U2);
	void compute_alternate_update_left(double* U1, double* U2);
	void compute_alternate_update_right(double* U1, double* U2);
	void compute_alternate_update_bottom(double* U1, double* U2);
	void compute_alternate_update_top(double* U1, double* U2);
	void compute_dist_squared();
	void compute_gen_sol();
	void cleanup();
};

#endif
