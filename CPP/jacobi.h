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
	void compute(int itermax,double e,bool neumannBC);
	void save();
	void save_gnuplot();
private:
// system variables
	operator_matrix *A;
	decomposition *D;
	comm_ctrl *C;
	double eps, dist_squared, sigma;
	double *RHS=NULL,*RHSit=NULL,*U=NULL,*Uit=NULL;
	int iterMax,iter, myNx,myNy,myN, Nx,Ny,N;
	bool NeumannBC;
// MPI variables
  int myRank, nOfProcs;
// file operations
	std::ofstream FILEOUT;
// functions
	void init(int itermax, double e, bool neumannBC);
	void init_sys();
	void init_MPI();
	void compute_iterate();
	void compute_iterate_local();
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
