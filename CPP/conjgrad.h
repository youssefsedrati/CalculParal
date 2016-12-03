#ifndef CONJGRAD_H
#define CONJGRAD_H

class CGmethod{
public:
	CGmethod();
	CGmethod(double Aii,double Cx,double Cy,int Nx,int N,double *RHS,double *U);
	~CGmethod();
	void compute(int iterMax, double eps);
private:
// system variables
	double Aii,Cx,Cy, eps,
					residue, drl, dwl, alpha, beta;
	double *RHS, *U, *r, *kappa, *d, *W;
	int iterMax,iter, Nx,Ny,N;
// MPI variables
  int myrank, nb_procs,
					fstline,lstline,fst,lst,deb,fin;
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
