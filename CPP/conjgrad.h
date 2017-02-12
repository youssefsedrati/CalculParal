#ifndef CONJGRAD_H
#define CONJGRAD_H

#include <iostream>
#include <fstream>
#include <string>
#include "operator_matrix.h"
#include "decomposition.h"
#include "comm_ctrl.h"

class CGMethod{
public:
	CGMethod();
	CGMethod(operator_matrix *a, decomposition *dc, double *RHS,double *U);
	~CGMethod();
	void compute(int iterMax, double eps);
	void save();
	void save_gnuplot();
private:
// system variables
	operator_matrix *A;
	decomposition *D;
	comm_ctrl *C;

	double eps, alpha, beta, gamma, norm;
	double *RHS=NULL,*RHSit=NULL,*U=NULL,*R=NULL,*P=NULL,*AP=NULL, *Z=NULL;
	int iterMax,iter, myNx,myNy,myN, Nx,Ny,N;
// MPI variables
  int myRank, nOfProcs;
// file operations
	std::ofstream FILEOUT;
// functions
	void init(int itermax, double e);
	void init_sys();
	void init_MPI();
	void compute_iterate();
	void compute_iterate_local();
	void compute_iterate_local_init();
	void compute_alpha();
	void compute_update();
	void compute_residue();
	void compute_beta();
	void compute_gradient();
	void compute_dist_squared();
	void compute_gen_sol();
	void cleanup();
	double vector_product(double *X, double *Y);
	double conjugate_vector_norm(double *X);
	void matrix_vector_product(double *X, double *RESULT);
	void matrix_vector_product_inner(double *X, double *RESULT);
	void matrix_vector_product_top(double *X, double *RESULT);
	void matrix_vector_product_bottom(double *X, double *RESULT);
	void matrix_vector_product_left(double *X, double *RESULT);
	void matrix_vector_product_right(double *X, double *RESULT);
	void matrix_vector_product_global(double *X, double *RESULT);
	void matrix_vector_product_global_inner(double *X, double *RESULT);
	void matrix_vector_product_global_top(double *X, double *RESULT);
	void matrix_vector_product_global_bottom(double *X, double *RESULT);
	void matrix_vector_product_global_left(double *X, double *RESULT);
	void matrix_vector_product_global_right(double *X, double *RESULT);
};


#endif
