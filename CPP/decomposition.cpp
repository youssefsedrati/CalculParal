#include "decomposition.h"

decomposition::decomposition(int myrank, int nb_procs, int nb_procs_x, int nx, int ny){
	myRank = myrank; N_procs = nb_procs; N_procs_x = nb_procs_x; N_procs_y = N_procs/N_procs_x;
	Nx = nx; Ny = ny;
	if(is_admissable){
		
	}
	else ~decomposition();
}

decomposition::~decomposition(){
	if(index_x) free(index_x);
	if(index_y) free(index_y);
}

bool decomposition::is_admissable(){
	return( (N_procs>0)&&(N_procs_x>0)&&(myRank<=N_procs)&&(N_procs_x*N_procs_y=N_procs) );
}