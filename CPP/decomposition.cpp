#include <stdlib.h>
#include "decomposition.h"

decomposition::decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny, bool iscovering = true){
	myRank = myrank; N_procs = nb_procs; N_procs_x = nb_procs_x; N_procs_y = N_procs/N_procs_x;
	Nx = nx; Ny = ny;
	isCovering = iscovering; // at this moment of time, true is the assumption, false case not yet handled
	if(is_admissable()){
		myRank_x = myRank%N_procs_x;
		myRank_y = (myRank-myRank_x)/N_procs_x;
		decompose_x();
		decompose_y();
		decompose_global();
	}
}

decomposition::~decomposition(){
	if(index_x) free(index_x);
	if(index_y) free(index_y);
}

bool decomposition::is_admissable(){
	return( (N_procs>0)&&(N_procs_x>0)&&(myRank<N_procs)&&(N_procs_x*N_procs_y==N_procs) );
}

void decomposition::decompose_x(){
	int begin = myRank_x*Nx/N_procs_x;
	if(N_procs_x==1) 
		myNx = Nx;
	else if(myRank_x==N_procs_x)
		myNx = Nx-begin;
	else
		myNx = 1+Nx/N_procs_x;
	index_x = (int*) malloc(myNx*sizeof(int));
	for(int i=0;i<myNx;++i){
		index_x[i] = begin+i;
	}
}

void decomposition::decompose_y(){
	index_y = (int*) malloc(myNy*sizeof(int));
}

void decomposition::decompose_global(){
	index_global = (int*) malloc(myN*sizeof(int));
}