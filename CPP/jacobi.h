#ifndef JACOBI_H
#define JACOBI_H

class JacobiMethod{
	JacobiMethod();
	JacobiMethod(double aii,double cx,double cy,int nx,int n,double *rhs,double *u);
	~JacobiMethod();
public:
	void compute(int itermax,double e);
private:
// system variables
	double Aii,Cx,Cy, eps, dist, sigma;
	double *RHS, *U,*Uit, *d, *W;
	int iterMax,iter, Nx,Ny,N;
// MPI variables
  int myrank, nb_procs,
					firstline,lastline,first,last,start,end;
  MPI_Status s1, s2;
  MPI_Request r1, r2, r3, r4;
// functions
	void init(int itermax, double e);
	void init_sys();
	void init_MPI();
	void compute_iterate();
	void compute_comm_main();
	void compute_alternate_update(double* U1, double* U2);
	void compute_gen_sol();
	void cleanup();
};

#endif
