#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "conjgrad.h"
#include "decomposition.h"
#include "operator_matrix.h"
#include "comm_ctrl.h"

/* implementation of the Conjugate Gradient.
	 iteration constants and vectors are successively updated in the routine
			compute_iterate();
	 each process uses its decomposition (ie. subdomain) to operate on;
*/

using namespace std;

// constructor & destructor
CGMethod::CGMethod(operator_matrix *a,decomposition *dc,
				double *rhs, double *u){
	A = a; D = dc; RHS = rhs; U = u;
	myN  = D->get_myN(); N  = D->get_N();
	myNx = D->get_myNx();Nx = D->get_Nx();
	myNy = D->get_myNy();Ny = D->get_Ny();
	Z  = (double*) malloc(myN*sizeof(double));
	R  = (double*) malloc(myN*sizeof(double));
  P  = (double*) malloc(myN*sizeof(double));
  AP  = (double*) malloc(myN*sizeof(double));
	RHSit= (double*) malloc(N*sizeof(double));
	C = new comm_ctrl(D,A,RHS,RHSit);
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
	//if(myRank!=0) return;
	cout << "#" << myRank << ". saving.\n";
	std::string filename = "CG_test_sol_";
	filename.append(to_string(myRank)); filename.append(".data");
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

void CGMethod::save_gnuplot(){ 
  if(myRank!=0) return; 
  cout << "#" << myRank << ". saving.\n"; 
  std::string filename = "conjgrad_gnuplot.dat"; 
  std::ofstream ofs("conjgrad_gnuplot.dat",std::ofstream::out); 
  if(!ofs) return; 
  for(int i=0;i<Ny;++i){ 
    for(int j=0;j<Nx;++j){ 
      double x= (double)(j+1)/(Nx+1), y= (double)(i+1)/(Ny+1); 
      ofs << x << " " << y << " " << U[j+Nx*i] << endl; 
    } 
    ofs << endl; 
  } 
  ofs.close(); 
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
	matrix_vector_product_global(U,AP);
	for(int i=0;i<N;++i){
		RHSit[i] = RHS[i];
		U[i] = 0;
	}
	int *idx = D->get_index_global(), j;
	for(int i=0;i<myN;++i){
		j = idx[i];
		U[j] = 1;
		R[i] = RHSit[j] - AP[i];
		P[i] = Z[i] = R[i]/A->Aii();
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
		C->send(U);// eventually wait until all processors have received
		C->cumulate_dist_squared(&norm);
		iter+=1;
	}
}

void CGMethod::compute_alpha(){
	gamma = vector_product(Z,R);
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
	matrix_vector_product_global(U,R);
	for(int i=0;i<myN;++i){
		int j = idx[i];
		R[i] = RHSit[j] - R[i];
		Z[i] = R[i]/A->Aii();
	}
	/*for(int i=0;i<myN;++i){
		R[i] = R[i] - alpha*AP[i];
	}*/
}

void CGMethod::compute_beta(){
	norm = vector_product(R,R);
	beta = norm/gamma;
}

void CGMethod::compute_gradient(){
	for(int i=0;i<myN;++i){
		P[i] = Z[i] + beta*P[i]/A->Aii();
	}
}

void CGMethod::compute_gen_sol(){
	MPI_Barrier(MPI_COMM_WORLD);
	C->compile_solution(U);
  if(myRank == 0){
    cout << "CG method terminated after " << iter-1
         << " iterations, squared_dist= "
         << sqrt(norm) << endl;
  }
}

void CGMethod::cleanup(){
	if(RHSit) free(RHSit);
	if(R) free(R);
	if(P) free(P);
	if(AP) free(AP);
}

double CGMethod::vector_product(double *X, double *Y){
	double sum = 0;
	for(int i=0;i<myN;++i)
		sum+= X[i]*Y[i];
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
	for(int j=1;j<myNy-1;++j){
		for(int k=1;k<myNx-1;++k){
			int i = k+j*myNx;
			RESULT[i]=A->Cx()*( X[i-1]+X[i+1] )
					 +A->Cy()*( X[i-myNx]+X[i+myNx] )
                     +A->Aii()* X[i];
		}
	}
}

void CGMethod::matrix_vector_product_top(double *X, double *RESULT){
	for(int j=1;j<myNx-1;++j){
		int i = j+(myNx)*(myNy-1);
		RESULT[i]=A->Cx()*( X[i-1]+X[i+1] )
				 +A->Cy()* X[i-myNx]
				 +A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_bottom(double *X, double *RESULT){
	for(int i=1;i<myNx-1;++i){
		RESULT[i]=A->Cx()*( X[i-1]+X[i+1] )
				 +A->Cy()* X[i+myNx]
				 +A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_left(double *X, double *RESULT){
	int i = 0;
	RESULT[i] = A->Cx()*X[i+1] + A->Cy()*X[i+myNx] + A->Aii()*X[i];
	i = myNx*(myNy-1);
	RESULT[i] = A->Cx()*X[i+1] + A->Cy()*X[i-myNx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		i = j*myNx;
		RESULT[i]=A->Cx()* X[i+1]
                 +A->Cy()*( X[i-myNx]+X[i+myNx] )
                 +A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_right(double *X, double *RESULT){
	int i = myNx-1;
	RESULT[i] = A->Cx()*X[i-1] + A->Cy()*X[i+myNx] + A->Aii()*X[i];
	i = myNx*myNy-1;
	RESULT[i] = A->Cx()*X[i-1] + A->Cy()*X[i-myNx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		int i = (j+1)*myNx-1;
		RESULT[i]=A->Cx()* X[i-1]
				 +A->Cy()*( X[i-myNx]+X[i+myNx] )
				 +A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_global(double *X, double *RESULT){
	matrix_vector_product_global_inner(X,RESULT);
	matrix_vector_product_global_top(X,RESULT);
	matrix_vector_product_global_bottom(X,RESULT);
	matrix_vector_product_global_left(X,RESULT);
	matrix_vector_product_global_right(X,RESULT);
}

void CGMethod::matrix_vector_product_global_inner(double *X, double *RESULT){
	int *idx=D->get_index_global_inner(), i;
	for(int j=1;j<myNy-1;++j){
		for(int k=1;k<myNx-1;++k){
			i = idx[(k-1)+(j-1)*(myNx-2)];
			RESULT[k+j*myNx] =A->Cx()*( X[i-1]+X[i+1] )
								+A->Cy()*( X[i-myNx]+X[i+myNx] )
								+A->Aii()* X[i];
		}
	}
}

void CGMethod::matrix_vector_product_global_top(double *X, double *RESULT){
	int *idx = D->get_index_global_top(), i;
	for(int j=1;j<myNx-1;++j){
		i = idx[j];
		RESULT[j+(myNx)*(myNy-1)] =A->Cx()*( X[i-1]+X[i+1] )
							+A->Cy()* X[i-myNx]
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_global_bottom(double *X, double *RESULT){
	int *idx = D->get_index_global_bottom(), i;
	for(int j=1;j<myNx-1;++j){
		i = idx[j];
		RESULT[j] =A->Cx()*( X[i-1]+X[i+1] )
							+A->Cy()* X[i+myNx]
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_global_left(double *X, double *RESULT){
	int *idx = D->get_index_global_left(), i;
	i = idx[0];
	RESULT[0] = A->Cx()*X[i+1] + A->Cy()*X[i+myNx] + A->Aii()*X[i];
	i = idx[myNy-1];
	RESULT[myNx*(myNy-1)] = A->Cx()*X[i+1] + A->Cy()*X[i-myNx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		i = idx[j];
		RESULT[j*myNx] =A->Cx()* X[i+1]
							+A->Cy()*( X[i-myNx]+X[i+myNx] )
							+A->Aii()* X[i];
	}
}

void CGMethod::matrix_vector_product_global_right(double *X, double *RESULT){
	int *idx = D->get_index_global_right(), i;
	i = idx[0];
	RESULT[myNx-1] = A->Cx()*X[i-1] + A->Cy()*X[i+myNx] + A->Aii()*X[i];
	i = idx[myNy-1];
	RESULT[myNx*myNy-1] = A->Cx()*X[i-1] + A->Cy()*X[i-myNx] + A->Aii()*X[i];
	for(int j=1;j<myNy-1;++j){
		i = idx[j];
		RESULT[(j+1)*myNx-1] =A->Cx()* X[i-1]
							+A->Cy()*( X[i-myNx]+X[i+myNx] )
							+A->Aii()* X[i];
	}
}

