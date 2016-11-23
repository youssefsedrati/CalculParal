#include <stdlib.h>
#include "decomposition.h"

decomposition::decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny){
	myRank = myrank; N_procs = nb_procs; N_procs_x = nb_procs_x; N_procs_y = N_procs/N_procs_x;
	Nx = nx; Ny = ny;
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
	if(index_global) free(index_global);
}

int* decomposition::get_index_x(){
	return index_x;
}

int* decomposition::get_index_y(){
	return index_y;
}

int* decomposition::get_index_global(){
	return index_global;
}

bool decomposition::is_admissable(){
	return( (N_procs>0)&&(N_procs_x>0)&&(myRank<N_procs)&&(N_procs_x*N_procs_y==N_procs) );
}

void decomposition::decompose_x(){
	int begin = myRank_x*Nx/N_procs_x;
	if(N_procs_x==1) 
		myNx = Nx;
	else if(myRank_x==N_procs_x)
		myNx = Nx-begin-1;
	else
		myNx = 1+Nx/N_procs_x;
	index_x = (int*) malloc(myNx*sizeof(int));
	for(int i=0;i<myNx;++i){
		index_x[i] = begin+i;
	}
}

void decomposition::decompose_y(){
	int begin = myRank_y*Ny/N_procs_y;
	if(N_procs_y==1) 
		myNy = Ny;
	else if(myRank_y==N_procs_y)
		myNy = Ny-begin-1;
	else
		myNy = 1+Ny/N_procs_y;
	index_y = (int*) malloc(myNy*sizeof(int));
	for(int i=0;i<myNy;++i){
		index_y[i] = begin+i;
	}
}

void decomposition::decompose_global(){
	myN = myNx * myNy;
	index_global = (int*) malloc(myN*sizeof(int));
	for(int i=0;i<myNx;++i)
		for(int j=0;j<myNy;++j)
			index_global[i+j*myNx]= index_x[i]+Nx*index_y[j];
}