#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "decomposition.h"

decomposition::decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny){
	myRank = myrank; N_procs = nb_procs; N_procs_x = nb_procs_x; N_procs_y = N_procs/N_procs_x;
	Nx = nx; Ny = ny; N = Nx*Ny;
	if(is_admissable()){
		myRank_x = myRank%N_procs_x;
		myRank_y = (myRank-myRank_x)/N_procs_x;
		decompose();
		accumulate_global_borders();
		accumulate_global_inner();
	}
	else std::cout << "Error: decomposition for #" << myrank << " is not admissable.\n";
}

decomposition::~decomposition(){
	if(index_x) free(index_x);
	if(index_y) free(index_y);
	if(index_global) free(index_global);
	if(index_global_top) free(index_global_top);
	if(index_global_bottom) free(index_global_bottom);
	if(index_global_left) free(index_global_left);
	if(index_global_right) free(index_global_right);
	if(index_global_inner) free(index_global_inner);
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

int* decomposition::get_index_global_top(){
	return index_global_top;
}

int* decomposition::get_index_global_bottom(){
	return index_global_bottom;
}

int* decomposition::get_index_global_left(){
	return index_global_left;
}

int* decomposition::get_index_global_right(){
	return index_global_right;
}

int* decomposition::get_index_global_inner(){
	return index_global_inner;
}

int decomposition::get_myNx(){
	return myNx;
}

int decomposition::get_myNy(){
	return myNy;
}

int decomposition::get_myN(){
	return myN;
}

int decomposition::get_myNinner(){
	return myNinner;
}

int decomposition::get_myRank_x(){
	return myRank_x;
}

int decomposition::get_myRank_y(){
	return myRank_y;
}
	
int decomposition::get_Nx(){
	return Nx;
}

int decomposition::get_Ny(){
	return Ny;
}

int decomposition::get_N(){
	return N;
}

int decomposition::get_N_procs(){
	return N_procs;
}

int decomposition::get_N_procs_x(){
	return N_procs_x;
}

int decomposition::get_N_procs_y(){
	return N_procs_y;
}

bool decomposition::is_admissable(){
	return( (N_procs>0)&&(N_procs_x>0)&&(myRank<N_procs)&&(N_procs_x*N_procs_y==N_procs) );
}

void decomposition::decompose(){
		decompose_x();
		decompose_y();
		decompose_global();
}

void decomposition::decompose_x(){
	int step = Nx/N_procs_x,
			carryover = Nx-N_procs_x*step,
			carryover_start = N_procs_x - carryover,
			begin = myRank_x*step;
	if(myRank_x==0) begin = 0;
	myNx = step;
	if(myRank_x+1>carryover_start){
		myNx++;
		begin = (carryover_start)*step+myNx*(myRank_x-carryover_start);
	}
	index_x = (int*) malloc(myNx*sizeof(int));
	for(int i=0;i<myNx;++i){
		index_x[i] = begin+i;
	}
}

void decomposition::decompose_y(){
	int step = Ny/N_procs_y,
			carryover = Ny-N_procs_y*step,
			carryover_start = N_procs_y - carryover,
			begin = myRank_y*step;
	if(myRank_y==0) begin = 0;
	myNy = step;
	if(myRank_y+1>carryover_start){
		myNy++;
		begin = (carryover_start)*step+myNy*(myRank_y-carryover_start);
	}
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

void decomposition::accumulate_global_borders(){
	accumulate_global_top();
	accumulate_global_bottom();
	accumulate_global_left();
	accumulate_global_right();
}

void decomposition::accumulate_global_top(){
	index_global_top = (int*) malloc(myNx*sizeof(int));
	int begin = (myNy-1)*myNx;
	for(int i=0;i<myNx;++i)
		index_global_top[i]= index_global[begin+i];
}

void decomposition::accumulate_global_bottom(){
	index_global_bottom = (int*) malloc(myNx*sizeof(int));
	for(int i=0;i<myNx;++i)
		index_global_bottom[i]= index_global[i];
}

void decomposition::accumulate_global_left(){
	index_global_left = (int*) malloc(myNy*sizeof(int));
	for(int i=0;i<myNy;++i)
		index_global_left[i]= index_global[i*myNx];		
}

void decomposition::accumulate_global_right(){
	index_global_right = (int*) malloc(myNy*sizeof(int));
	for(int i=0;i<myNy;++i)
		index_global_right[i]= index_global[(i+1)*myNx-1];
}

void decomposition::accumulate_global_inner(){
	myNinner = (myNx-2)*(myNy-2);
	if(myNinner<1) return;
	index_global_inner = (int*) malloc(myNinner*sizeof(int));
	for(int i=1;i<myNx-1;++i)
		for(int j=1;j<myNy-1;++j)
			index_global_inner[i-1+(j-1)*(myNx-2)]= index_global[i+j*myNx];
}