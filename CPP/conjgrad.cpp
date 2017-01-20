#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "conjgrad.h"
#include "decomposition.h"
#include "operator_matrix.h"
#include "comm_ctrl.h"

using namespace std;

// constructor & destructor
CGMethod::CGMethod(operator_matrix *a,decomposition *dc,
				double *rhs, double *u){
	A = a; D = dc; RHS = rhs; U = u;
	myN  = D->get_myN(); N  = D->get_N();
	myNx = D->get_myNx();Nx = D->get_Nx();
	myNy = D->get_myNy();Ny = D->get_Ny();
	Uup  = (double*) malloc(N*sizeof(double));
  R  = (double*) malloc(myN*sizeof(double));
  P  = (double*) malloc(myN*sizeof(double));
  AP  = (double*) malloc(myN*sizeof(double));
	RHSit= (double*) malloc(N*sizeof(double));
	C = new comm_ctrl(D,A,RHS,RHSit,Uup);
}

CGMethod::~CGMethod(){
	cleanup();
	delete C;
}
// public
void CGMethod::compute(int itermax, double e){
	init(itermax,e);
	compute_iterate();	
	compute_gen_sol();
}

void CGMethod::save(){
	if(myRank!=0) return;
	cout << "#" << myRank << ". saving.\n";
	std::string filename = "Jacobi_test_sol.data";
	//filename.append(to_string(myRank)); filename.append(".data");
	FILEOUT = ofstream(filename,ios::out);
	if(!FILEOUT) return;
	FILEOUT << Nx << " "<< Ny <<endl;
	for(int i=Ny-1;i>=0;--i){
		for(int j=0;j<Nx;++j)
			FILEOUT << U[j+Nx*i] << " ";
		FILEOUT << endl;
	}
	FILEOUT.close();
}

// private
void CGMethod::init(int itermax, double e){
	iterMax = itermax; eps = e;
  iter = norm = 1;
	init_MPI();
	init_sys();
}

void CGMethod::init_MPI(){	
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nOfProcs);
}

void CGMethod::init_sys(){
	matrix_vector_product(U,AP);
	for(int i=0;i<N;++i){	
		RHSit[i] = RHS[i];
		U[i] = 0;
	}
	int *idx = D->get_index_global(), j;
	for(int i=0;i<myN;++i){
		j = idx[i];
		P[j] = R[j] = RHS[j] - AP[j];
	}
}

void CGMethod::compute_iterate(){	
	while( (iter<=iterMax)&&(norm>eps*eps) ){
		C->receive(); // eventually wait until all processors have sent
		compute_alpha();
		compute_update();
		compute_residue();
		compute_beta();
		compute_gradient();
		C->send(U);// eventually wait until all processors have receivedte_dist_squared();
		C->cumulate_dist_squared(&norm);
		iter+=1;
	}
}

void CGMethod::compute_alpha(){
	gamma = vector_product(R,R);
	alpha = gamma/conjugate_vector_norm(P);
}

void CGMethod::compute_update(){
	int *idx = D->get_index_global();
	for(int j=0;j<myN;++j){
		int i=idx[j];
		U[i] = U[i] + alpha*P[j];
	}
}

void CGMethod::compute_residue(){
	int *idx = D->get_index_global();
	for(int j=0;j<myN;++j){
		int i=idx[j];
		R[i] = R[i] - alpha*AP[i];
	}
}

void CGMethod::compute_beta(){
	norm = vector_product(R,R);
	beta = norm/gamma;
}

void CGMethod::compute_gradient(){
	int *idx = D->get_index_global();
	for(int j=0;j<myN;++j){
		int i=idx[j];
		P[i] = R[i] + beta*P[i];
	}
}

void CGMethod::compute_gen_sol(){
	MPI_Barrier(MPI_COMM_WORLD);
	C->compile_solution(U);		
  if(myRank == 0){
    printf("CG method terminated after"
           " %d iterations, squared_dist= %0.31f\n", iter-1,sqrt(norm));
  }
}

void CGMethod::cleanup(){
	if(RHSit) free(RHSit);
	if(R) free(R);
	if(P) free(P);
	if(AP) free(AP);
	if(Uup) free(Uup);
}

double CGMethod::vector_product(double *X, double *Y){
	double sum = 0;
	int *idx = D->get_index_global();
	for(int j=0;j<myN;++j){
		int i=idx[j];
		sum+= X[i]*Y[i];
	}
	return sum;
}

double CGMethod::conjugate_vector_norm(double *X){
	matrix_vector_product(X,AP);
	return vector_product(X,AP);
}

void CGMethod::matrix_vector_product(double *X, double *RESULT){
	matrix_vector_product_inner(X,RESULT);
	matrix_vector_product_top(X,RESULT);
	matrix_vector_product_bottom(X,RESULT);
	matrix_vector_product_left(X,RESULT);
	matrix_vector_product_right(X,RESULT);
}

void CGMethod::matrix_vector_product_inner(double *X, double *RESULT){
	int *idx = D->get_index_global_inner();
	for(int j=0;j<D->get_myNinner();++j){
		int i = idx[j];
		RESULT[i] =A->Cx()*( X[i-1]+X[i+1] )
							+A->Cy()*( X[i-Nx]+X[i+Nx] )
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_top(double *X, double *RESULT){
	int *idx = D->get_index_global_top();
	for(int j=1;j<myNx-1;++j){
		int i = idx[j];
		RESULT[i] =A->Cx()*( X[i-1]+X[i+1] )
							+A->Cy()* X[i-Nx]
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_bottom(double *X, double *RESULT){
	int *idx = D->get_index_global_bottom();
	for(int j=1;j<myNx-1;++j){
		int i = idx[j];
		RESULT[i] =A->Cx()*( X[i-1]+X[i+1] )
							+A->Cy()* X[i+Nx]
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_left(double *X, double *RESULT){
	int *idx = D->get_index_global_left(),
	i = idx[0];
	RESULT[i] = A->Cx()*X[i+1] + A->Cy()*X[i+Nx] + A->Aii()*X[i];
	i = idx[myNy-1];
	RESULT[i] = A->Cx()*X[i+1] + A->Cy()*X[i-Nx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		i = idx[j];
		RESULT[i] =A->Cx()* X[i+1] 
							+A->Cy()*( X[i-Nx]+X[i+Nx] )
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_right(double *X, double *RESULT){
	int *idx = D->get_index_global_right(),
	i = idx[0];
	RESULT[i] = A->Cx()*X[i-1] + A->Cy()*X[i+Nx] + A->Aii()*X[i];
	i = myNy-1;
	RESULT[i] = A->Cx()*X[i-1] + A->Cy()*X[i-Nx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		int i = j*myNx-1;
		RESULT[i] =A->Cx()* X[i-1] 
							+A->Cy()*( X[i-Nx]+X[i+Nx] )
							+A->Aii()* X[i];
	}
}